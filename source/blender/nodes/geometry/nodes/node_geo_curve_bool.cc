/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "BKE_curves.hh"
#include "BKE_spline.hh"
#include "BLI_timeit.hh"
#include "UI_interface.h"
#include "UI_resources.h"
#include "node_geometry_util.hh"
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <list>

using blender::attribute_math::mix2;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

/**
 * TODO :
 * Turn Curves into CurveHull asap.
 * Then work just with CurveHull.
 * CurveHull should never need to refer to the actual curve again after creation.
 *
 */

/**
 * This file implements the Curve 2D Boolean geometry node.
 * The recurring variable v exclusively describes a [0,1] float that represents the distance along
 * a curve or a line. All other names should be self explanatory.
 *
 *
 * EXPLANATION :
 *
 * BOOLEAN ALGORITHM :
 *
 * (insert video url here)
 * In this video I'm going to explain to you how the Boolean Curve Algorithm I wrote for blender
 * works.
 * ( Show Simple overlapping circles, with A and B inside )
 * For the time being, let's assume all curves go clockwise, curves don't intersect themselves, and
 * all curves are closed. And let's try doing a UNION operation first.
 * Each curve exists as a doubly-linked-list of points. Meaning
 * we can iterate them in linear time, insert and remove in constant time, and find in linear time.
 * Let's start by calculating all intersection points and inserting them into both curve A and
 * curve B. Meaning each intersection point generates 2 new curve points, one on each curve, but
 * each at the same position. These intersection points will point to one another, turning the
 * doubly-linked-list of each curve into a partial tripply-linked-mesh of all intersecting curves.
 * ( Show Intersection points )
 * Assuming we found intersection points, let's take the first point on curve A.
 * There are 2 cases here.
 * 1) The point lies outside curve B.
 * 2) The point lies inside curve B.
 * We then trace along the curve until we find a point that is an intersection point.
 * Now let's assume we're doing a UNION operation. In that case, any point on A past the
 * intersection point will be outside of the other shape if it was inside before, or inside if it
 * was outside before. Since UNION just means we want our final shape to take up the maximum area
 * it can, we want to follow whatever curve is outside. In this case, A was inside before, so we
 * discard anything we traced before, and continue tracing on curve A, because A is now OUTSIDE. We
 * then reach the second intersection point. At this point curve A was outside before, meaning
 * after this point it will be inside. So we switch to curve B, which must now be outside. We
 * continue this until we reach a point we've already processed. We then repeat the algorithm with
 * any curve points we did not process.
 *
 * Things go very similarly for the other 2 operations; INTERSECTION and DIFFERENCE.
 * For DIFFERENCE, we just need to flip the direction of our subtracting curve. This means that
 * when we flip from curve A to curve B, we move INSIDE curve A along curve B instead of OUTSIDE.
 * For INTERSECTION we just flip our "is_inside" function results, since we want to MINIMIZE area
 * now. This "is_inside" value is also great to form a mathematical proof on, but I'm not gonna go
 * into that here. I recommend just taking a piece of paper and trying this by hand. It makes
 * intuitive sense.
 *
 * From the DIFFERENCE operation we've now learned that a subtractive shape is nothing but a union
 * with a curve that goes the opposite way it's supposed to. This is a generalizable property. We
 * know that the shape of a donut ,for example, has an inside and an outside curve. The outside
 * curve describes where the shape starts, while the inside curve describes an empty area within
 * the shape. Putting these 2 concepts together, we can process complicated shapes with multiple
 * curves by flipping the direction of inner curves. Every curve that is surrounded by an odd
 * number of curves is an inner (subtractive) curve.
 *
 *
 * 2D CONCAVE INTERSECTION
 *
 * TODO: What runtime is actually feasable here? What about prep time?
 * The goal of this algorithm is to find out if 2 concave polygons intersect.
 * The target runtime is O( (log(n) + log(m))^2 ) where n and m are the number of points in the
 * respective polygon. The target prep time is O( n * log(n) ) for each polygon. We do this by
 * first splitting the concave polygons into n convex parts. ( O(n) ) We then build a convex hull
 * binary tree. O( n * log(n) )
 *
 * Worst case : There are n/2 and m/2 convex segments that all overlap. Runtime would be O( n * m
 * ). So let's look at worst case with NO overlaps instead.
 *
 *
 * 2D CONVEX INTERSECTION
 *
 * The goal of this algorithm is to find out if 2 convex polygons intersect.
 * The target runtime is O( log(n) + log(m) ) where n and m are the number of points in the
 * respective polygon.
 *
 *
 *
 * BUGS:
 *
 * TODO:
 * Turn as many stds into blender collections as possible.
 * Move static intersection to BezierSpline
 *
 * Optimize multiple bool operations in a row via multiple operands?
 * Line Segment intersection library
 * Optional hard corners so attributes like radius stick. Or maybe an enum to define which curve
 * should provide attributes for intersection points. Conformity slider : If 1, map intersection
 * points to geometry, if 0, map to bezier curve. Interpolate between. If either of the curves self
 * intersect, show a warning.
 * Splines to bezier conversion. ( How does this work with weighted points? We need 1 for 1 bezier
 * control points, no?)
 */

/**
 * Same ids as "curve_bool_items" in NOD_static_types.h
 */
typedef enum BoolOperationType {
  OR = 0,
  AND = 1,
  SUB = 2,
} BoolOperationType;

/**
 * DEBUG CODE START
 */

static void print(std::string s)
{
  std::cout << s << std::endl << std::flush;
}

template<typename T> static void print(T s)
{
  std::cout << s << std::endl << std::flush;
}

template<typename T> static std::string str(T t)
{
  return std::to_string(t);
}

/**
 * DEBUG CODE END
 */

namespace blender::nodes {

/**
 * This is called whenever the node needs to update, such as when inputs change.
 */
static void geo_node_curve_bool_exec(GeoNodeExecParams params)
{
  SCOPED_TIMER("Curve Boolean");

  int resolution = params.extract_input<int>("Int. Resolution");
  GeometrySet geometry_set_a = params.extract_input<GeometrySet>("Base Curve");
  // CurveComponent& curve_set_a =
  // geometry_set_a.get_component_for_write<CurveComponent>();//bke::geometry_set_realize_instances(curve_set_a);
  std::vector<const CurveEval *> curves_b;

  Vector<GeometrySet> geometry_sets = params.extract_multi_input<GeometrySet>("Curves");

  if (!(geometry_set_a.has_curves())) {
    print("Nope");
    params.set_output("Curve", geometry_set_a);
    return;
  }

  const bNode &node = params.node();
  const BoolOperationType data_type = static_cast<BoolOperationType>(node.custom2);

  int real_count = 0;
  // Repeat the bool operation for each element in the operand `Curves` input.
  Curves *primary_curves = geometry_set_a.get_curves_for_write();
  std::unique_ptr<Curves> result_curves;
  for (GeometrySet &set_group : geometry_sets) {
    // const CurveComponent* curve_in = set_group.get_component_for_read<CurveComponent>();
    if (set_group.has_curves()) {
      // result_curves = generate_boolean_shapes(
      //     primary_curves, set_group.get_curves_for_write(), data_type, resolution);

      primary_curves = result_curves.get();

      real_count++;
    }
  }

  if (real_count == 0) {
    params.set_output("Curve", geometry_set_a);
  }
  else {
    params.set_output("Curve", GeometrySet::create_with_curves(result_curves.release()));
  }
}

static void geo_node_curve_bool_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Geometry>("Base Curve");
  b.add_input<decl::Geometry>("Curves").multi_input();
  b.add_input<decl::Int>("Int. Resolution")
      .min(1)
      .description(
          "If multiple bezier curves are used, this node will process them one by one and "
          "generate an intermediate curve after each step. Lowering that curve's resolution will "
          "increase performance but decrease accuracy. This does not affect input curves. It does "
          "affect the output.")
      .default_value(12);
  b.add_output<decl::Geometry>("Curve");
}

/**
 * This function is needed to draw the enum dropdown of operations.
 */
static void geo_node_curve_bool_layout(uiLayout *layout, bContext *UNUSED(C), PointerRNA *ptr)
{
  uiLayoutSetPropSep(layout, true);
  uiLayoutSetPropDecorate(layout, false);
  uiItemR(layout, ptr, "operation", 0, "", ICON_NONE);
}

}  // namespace blender::nodes

void register_node_type_geo_curve_bool()
{
  static bNodeType ntype;
  geo_node_type_base(&ntype, GEO_NODE_CURVE_BOOL, "Curve Bool", NODE_CLASS_GEOMETRY);
  ntype.geometry_node_execute = blender::nodes::geo_node_curve_bool_exec;
  ntype.draw_buttons = blender::nodes::geo_node_curve_bool_layout;
  ntype.declare = blender::nodes::geo_node_curve_bool_declare;
  nodeRegisterType(&ntype);
}
