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
 * - clone original data.
 * - get curves as indices/offsets list.
 * - calculate intersections. Use list iterator, so insert() can add it right after. Add data to
 * end of cloned data, with full interpolation.
 * - trace.
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
 * \warning Copied from subdivide_curves.cc
 * In theory this should copy the curves, but all I get is an exception down the line.
 */
static Curves *create_result_curves(const bke::CurvesGeometry &src_curves)
{
  Curves *dst_curves_id = bke::curves_new_nomain(0, src_curves.curves_num());
  bke::CurvesGeometry &dst_curves = bke::CurvesGeometry::wrap(dst_curves_id->geometry);
  CurveComponent dst_component;
  dst_component.replace(dst_curves_id, GeometryOwnershipType::Editable);
  /* Directly copy curve attributes, since they stay the same. */
  CustomData_copy(&src_curves.curve_data,
                  &dst_curves.curve_data,
                  CD_MASK_ALL,
                  CD_DUPLICATE,
                  src_curves.curves_num());
  dst_curves.runtime->type_counts = src_curves.runtime->type_counts;

  return dst_curves_id;
}

/**
 * Generate one or more new curves from 2 existing sets of curves.
 * These curves must not self intersect.
 * The general idea is to follow one of the curves and copy the control points until an
 * intersection is found. At which point we know the other curve is bigger, so we switch to the
 * other one. We do this for each intersection until we reach our initial position.
 */
static Curves* generate_boolean_shapes(GeometrySet &primary_geometry_set,
                                                       GeometrySet &secondary_geometry_set,
                                                       BoolOperationType type,
                                                       int resolution)
{

  const CurveComponent &primary_component =
      *primary_geometry_set.get_component_for_read<CurveComponent>();
  const Curves &_primary_curves = *primary_component.get_for_read();
  bke::CurvesGeometry primary_curves = bke::CurvesGeometry::wrap(_primary_curves.geometry);
  const int primary_curve_count = primary_curves.curve_num;

  // const CurveComponent &secondary_component =
  //     *secondary_geometry_set.get_component_for_read<CurveComponent>();
  // const Curves &_secondary_curves = *secondary_component.get_for_read();
  // bke::CurvesGeometry secondary_curves = bke::CurvesGeometry::wrap(_secondary_curves.geometry);
  // const int secondary_curve_count = secondary_curves.curve_num;

  // std::vector<int> results = {};

  // // TODO : Walk along all curves in A, put the offsets in array. Then turn array into valid output
  // for (const int i_curve : primary_curves.curves_range()) {
  //   const IndexRange curve_index_range =  primary_curves.points_for_curve(i_curve);
  //   for (const int i_point : curve_index_range) {
  //     results.push_back(i_point);
  //   }
  // }

  return create_result_curves(primary_curves);
}

/**
 * This is called whenever the node needs to update, such as when inputs change.
 */
static void geo_node_curve_bool_exec(GeoNodeExecParams params)
{
  SCOPED_TIMER("Curve Boolean");

  int resolution = params.extract_input<int>("Int. Resolution");

  GeometrySet geometry_set_a = params.extract_input<GeometrySet>("Base Curve");

  Vector<GeometrySet> geometry_sets = params.extract_multi_input<GeometrySet>("Curves");

  if (!(geometry_set_a.has_curves())) {
    print("Nope");
    params.set_output("Curve", geometry_set_a);
    return;
  }

  const bNode &node = params.node();
  const BoolOperationType data_type = static_cast<BoolOperationType>(node.custom2);

  int real_count = 0;
  // Repeat the bool operation for each element in the `Curves` input.
  Curves *primary_curves = geometry_set_a.get_curves_for_write();
  std::unique_ptr<Curves> result_curves;
  for (GeometrySet &set_group : geometry_sets) {
    // const CurveComponent* curve_in = set_group.get_component_for_read<CurveComponent>();
    if (set_group.has_curves()) {
      primary_curves = generate_boolean_shapes(geometry_set_a, set_group, data_type, resolution);

      real_count++;
    }
  }

  if (real_count == 0) {
    params.set_output("Curve", geometry_set_a);
  }
  else {
    params.set_output("Curve", GeometrySet::create_with_curves(primary_curves));
  }

}

/**
 * Define the shape of the node, its inputs and outputs.
 */
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
