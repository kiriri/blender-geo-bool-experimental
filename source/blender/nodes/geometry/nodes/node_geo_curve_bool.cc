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
using HandleType = BezierSpline::HandleType;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

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

// TODO : There's got to be a blender version of this...
struct Rect {
  blender::float2 min;
  blender::float2 max;

  Rect(blender::float2 min, blender::float2 max) : min(min), max(max)
  {
  }

  /**
   * Check if the other #Rect intersects with this one.
   */
  bool intersect(Rect other)
  {
    return (this->min.x < other.max.x) && (this->max.x > other.min.x) &&
           (this->min.y < other.max.y) && (this->max.y > other.min.y);
  }

  /**
   * Same as intersect, but componentwise.
   * Returns 0 if no overlap occured, 1 if x axis overlapped, 2 if y axis overlapped, 3 if true
   * intersection occured
   */
  bool partial_intersect(Rect other)
  {
    return ((this->min.x < other.max.x) && (this->max.x > other.min.x)) +
           (((this->min.y < other.max.y) && (this->max.y > other.min.y)) << 1);
  }

  /**
   * Check if a point lies within this rect.
   */
  bool contains(blender::float2 point)
  {
    return (point.x <= this->max.x) && (point.x >= this->min.x) && (point.y <= this->max.y) &&
           (point.y >= this->min.y);
  }

  /**
   * Check if a #Rect is completely enveloped by this one.
   */
  bool contains(Rect other)
  {
    return (this->min.x < other.min.x) && (this->max.x > other.max.x) &&
           (this->min.y < other.min.y) && (this->max.y > other.max.y);
  }

  /**
   * Expand this to contain the other #Rect.
   */
  void expand(Rect other)
  {
    this->min.x = min_ff(this->min.x, other.min.x);
    this->max.x = max_ff(this->max.x, other.max.x);
    this->min.y = min_ff(this->min.y, other.min.y);
    this->max.y = max_ff(this->max.y, other.max.y);
  }

  /**
   * Create a new #Rect that contains both this and the other #Rect.
   */
  Rect Combined(Rect other)
  {
    return {blender::float2(min_ff(this->min.x, other.min.x), min_ff(this->min.y, other.min.y)),
            blender::float2(max_ff(this->max.x, other.max.x), max_ff(this->max.y, other.max.y))};
  }

  float area()
  {
    float x = this->max.x - this->min.x;
    float y = this->max.y - this->min.y;
    return sqrtf(x * x + y * y);
  }

  float area_sqr()
  {
    float x = this->max.x - this->min.x;
    float y = this->max.y - this->min.y;
    return x * x + y * y;
  }

  bg::model::box<bg::model::point<float, 2, bg::cs::cartesian>> Box()
  {
    bg::model::point<float, 2, bg::cs::cartesian> min = {this->min.x, this->min.y};
    bg::model::point<float, 2, bg::cs::cartesian> max = {this->max.x, this->max.y};
    return bg::model::box<bg::model::point<float, 2, bg::cs::cartesian>>(min, max);
  }
};

// namespace blender::trees {

// static const int NODE_CAPACITY = 3;
// template<typename T> struct RTree_Node {
//   Rect rect;
//   T value;
//   RTree_Node *children[NODE_CAPACITY];
//   // How many children are in use.
//   int size = 0;
//   // How many elements in total exist in this subtree.
//   int real_size = 1;

//   /**
//    * Try to move any of the current children under another one if the new rect is more suitable
//    to
//    * exist as a child of this node.
//    * Examples of this would be if one of the current children can be fully or mostly contained
//    * within another current child.
//    * Returns whether space could be made. In which case r_index will return which index was
//    cleared.
//    * Otherwise r_index will return the most suitable index for the rect to be inserted under.
//    */
//   bool _make_space(Rect rect, int &r_index)
//   {
//     // If there's no space for another leaf, either one of the existing subtrees or the new
//     value
//     // needs to be put in a subtree
//     float rect_area = rect.area_sqr();
//     float min_wasted_area = INFINITY;
//     RTree_Node *min_child = nullptr;

//     for (int i = 0; i < size; i++) {
//       Rect child_rect = this->children[i]->rect;

//       Rect combined_rect = child_rect.Combined(rect);
//       float combined_area = combined_rect.area_sqr();
//       float child_area = child_rect.area_sqr();
//       // this is inaccurate, because overlap will count twice (which is sorta good in this
//       // heuristic)
//       float wasted_area = combined_area - rect_area - child_area;

//       //  if any of the existing rects contains the new one entirely, put it as a child.
//       if (child_rect.contains(rect)) {

//         this->children[i]->insert(rect, value);
//         return;
//       }

//       if (min_wasted_area < wasted_area) {
//         min_wasted_area = wasted_area;
//         min_child = this->children[i];
//       }
//     }

//     // If it can be a direct child and it doesn't fit with any of the existing children
//     if (this->size < NODE_CAPACITY && min_wasted_area > rect_area) {
//     }
//   }

//   /**
//    * Insert a value under a specific #Rect area.
//    */
//   void insert(Rect rect, T value)
//   {
//     this->real_size++;

//     // If this node has the size, insert it. Even if it is theoretically contained by another
//     // child.
//     if (this->size < NODE_CAPACITY) {
//       // If this is a leaf, turn it into a branch and create the original as a leaf.
//       if (this->size == 0) {
//         this->children[size++] = {this->rect, this->value};
//       }
//       this->children[size++] = {rect, value};
//       this->rect.expand(rect);
//       return;
//     }

//     this->rect.expand(rect);

//     // If any of the
//     bool contained_one = false;
//     for (int i = 0; i < size; i++) {
//       Rect child_rect = this->children[i]->rect;
//       if (rect.contains(child_rect)) {

//         // Keep going to check if more are contained.
//         contained_one = true;
//       }
//     }
//     if(contained_one)
//       return;

//     int index;
//     this->_make_space(rect, index);
//   }

//   /**
//    * Returns all #Rect elements that intersect with
//    */
//   void find(std::vector<Rect> r_results)
//   {
//   }
// };

// /**
//  * Custom 2D R-Tree implementation.
//  * R-Trees are just hierarchical bounding boxes.
//  * They are not balanced.
//  * Worst case :
//  * In the case where for all Rects N(x) in N, N(x-1) contains N(x), this produces a linear list.
//  */
// template<typename T> class RTree {
//   RTree_Node<T> *root;

//   // TODO : Could we theoretically keep this balanced?
//   // If so, we could allocate a vector

//  public:
//   void insert(Rect rect, T value)
//   {
//     if (root == nullptr) {
//       root = new RTree_Node<T>{rect, value};
//     }
//     else {
//       root->insert(rect, value);
//     }
//   }
// };
// }  // namespace blender::trees

namespace blender::nodes {
/**
 * Does the curve contain the point?
 * If the curve is open, it will be interpreted as a closed shape with a straight line between
 * start and end.
 */
static bool curve_contains(Spline *spline, float3 point)
{
  // Using bounds would not improve performance, because they aren't cached.
  blender::Span<float3> points = spline->evaluated_positions();
  int points_length = points.size();
  int i, j;
  bool result = 0;
  for (i = 0, j = points_length - 1; i < points_length; j = i++) {
    if (((points[i].y > point.y) != (points[j].y > point.y)) &&
        (point.x <
         (points[j].x - points[i].x) * (point.y - points[i].y) / (points[j].y - points[i].y) +
             points[i].x)) {
      result = !result;
    }
  }

  return result;  // if even, point is outside ( return false )
}

/**
 * Count the number of splines that contain the given point.
 * Ignore up to one particular spline, usually the one the point is located on
 */
static int count_containing_curves(Span<Spline *> splines,
                                   float3 point,
                                   const Spline *except = nullptr)
{
  int counter = 0;
  for (Spline *spline : splines) {
    if ((spline != except) && curve_contains(spline, point)) {
      counter++;
    }
  }
  return counter;
}

/**
 * Find intersection between 2 2d line segments, each defined by a start and end point. Z
 * Coordinates are ignored. If either of the 2 floats returned is negative, no collision occured.
 * If they are not negative, they are in [0,1]
 */
static float2 linesIntersect(float2 A, float2 B, float2 C, float2 D)
{
  float AC_x = C.x - A.x;
  float AC_y = C.y - A.y;
  float AB_x = B.x - A.x;
  float AB_y = B.y - A.y;
  float CD_x = D.x - C.x;
  float CD_y = D.y - C.y;

  float ACxAB = AC_x * AB_y - AC_y * AB_x;

  // Lines are collinear. They intersect infinitely if they have any overlap. We can't handle that.
  if (ACxAB == 0.f) {
    return float2(-1, -1);
  }

  float ACxCD = AC_x * CD_y - AC_y * CD_x;
  float ABxCD = AB_x * CD_y - AB_y * CD_x;
  // Multiplication is MUCH faster than division. Worth the extra operation.
  float one_over_ABxCD = 1.f / ABxCD;
  float t = ACxCD * one_over_ABxCD;
  float u = ACxAB * one_over_ABxCD;

  if ((t >= 0.f) && (t <= 1.f) && (u >= 0.f) && (u <= 1.f)) {
    return float2(t, u);
  }
  return float2(-1, -1);
}

static bool linesIntersect_b(float2 A, float2 B, float2 C, float2 D)
{
  float2 results = linesIntersect(A, B, C, D);
  return results.x >= 0 && results.y >= 0;
}

/**
 * Same as linesIntersect, but with rays.
 * Input must not be parallel.
 */
static float2 raysIntersect(float2 A, float2 B, float2 C, float2 D)
{
  float AC_x = C.x - A.x;
  float AC_y = C.y - A.y;
  float AB_x = B.x - A.x;
  float AB_y = B.y - A.y;
  float CD_x = D.x - C.x;
  float CD_y = D.y - C.y;

  float ACxAB = AC_x * AB_y - AC_y * AB_x;
  float ACxCD = AC_x * CD_y - AC_y * CD_x;
  float ABxCD = AB_x * CD_y - AB_y * CD_x;
  // Multiplication is MUCH faster than division. Worth the extra operation.
  float one_over_ABxCD = 1.f / ABxCD;
  float t = ACxCD * one_over_ABxCD;
  float u = ACxAB * one_over_ABxCD;

  return float2(t, u);
}

static float2 raysIntersectSimple(float2 A, float2 B, float2 C, float2 D)
{
  float AC_x = C.x - A.x;
  float AC_y = C.y - A.y;
  float AB_x = B.x - A.x;
  float AB_y = B.y - A.y;
  float CD_x = D.x - C.x;
  float CD_y = D.y - C.y;

  float ACxCD = AC_x * CD_y - AC_y * CD_x;
  float ABxCD = AB_x * CD_y - AB_y * CD_x;

  float t = ACxCD / ABxCD;

  return float2(A.x + AB_x * t, A.y + AB_y * t);
}

struct CurveHull;

struct Intersection {
  // distance between the current evaluated point and the next, in [0,1)
  float v;
  // Intersection on the opposite segment
  Intersection *target;
  // self
  CurveHull *hull;
};

/**
 * Stores a complete map of all connected points and intersections across all involved splines as a
 * double linked list.
 */
struct CurveHull {
  // Used to identify which curve we're on, which is useful if self-intersection is disabled.
  int curve_id;
  // original spline this point is on.
  Spline *spline;
  // Index of evaluated point on curve.
  int curve_i;
  // Distance along curve, in [0,1) .
  float v;
  // Real evaluated position on the original curve.
  float3 position;
  // Next point along the original curve.
  CurveHull *next = nullptr;
  // Previous point along the original curve.
  CurveHull *prev = nullptr;
  // If not nullptr, this is a virtual CurveHull point, that denotes an intersection.
  Intersection *intersection = nullptr;
  // Used to collect intersections without creating additional segments. Only control points have
  // anything in here.
  std::vector<Intersection *> intersections;
  // Used when tracing to determine if a point has been passed through.
  int pass = -1;
  // Used when tracing. If true, the next intersection is passed right through.
  // 0 = false, 1 = true, -1 = undefined. Only starting points have this defined.
  int is_inside = -1;

  /**
   * Find the next control point that HAS intersections.
   */
  CurveHull *next_intersection_control(bool include_self = false, CurveHull *__start = nullptr)
  {
    if (include_self && ((intersections.size() > 0))) {
      return this;
    }

    if (next == nullptr || __start == this) {  // end condition
      return nullptr;
    }

    if (__start == nullptr) {
      __start = this;
    }

    return next->next_intersection_control(true, __start);
  }

  /**
   * Find the next point that IS an intersection.
   */
  CurveHull *next_intersection(bool include_self = false, CurveHull *__start = nullptr)
  {
    if (include_self && intersection != nullptr) {
      return this;
    }

    if (next == nullptr || __start == this) {  // end condition
      return nullptr;
    }

    if (__start == nullptr) {
      __start = this;
    }

    return next->next_intersection(true, __start);
  }

  /**
   * Removes this and all following points as well as their intersection structs.
   */
  void remove()
  {
    if (prev != nullptr) {
      prev->next = nullptr;
    }

    if (next != nullptr) {
      next->prev = nullptr;
      next->remove();
    }

    if (intersection != nullptr) {
      delete intersection;
    }

    delete this;
  }
};

/**
 * Convert a path into a CurveHull form, which is a double linked list with extra variables for
 * tracing along intersections lateron.
 */
static CurveHull *points_to_hull(Spline *spline, int curve_id)
{
  blender::Span<float3> points = spline->evaluated_positions();
  bool cycle = spline->is_cyclic();

  CurveHull *start = new CurveHull{.curve_id = curve_id,
                                   .spline = spline,
                                   .curve_i = 0,
                                   .v = 0,
                                   .position = points[0],
                                   .next = nullptr,
                                   .prev = nullptr};
  CurveHull *prev = start;

  int size = points.size();
  int real_size = size + (cycle ? 0 : -1);
  for (int i = 1; i < size; i++) {
    CurveHull *current = new CurveHull{.curve_id = curve_id,
                                       .spline = spline,
                                       .curve_i = i,
                                       .v = i / (float)real_size,
                                       .position = points[i],
                                       .prev = prev};
    prev->next = current;
    prev = current;
  }

  if (cycle) {
    prev->next = start;
    start->prev = prev;
  }

  return start;
}

// static bool any_collide(std::vector<ConvexHull *> convex_hulls_a,
//                         std::vector<ConvexHull *> convex_hulls_b)
// {
//   for (auto convex_a : convex_hulls_a) {
//     for (auto convex_b : convex_hulls_b) {
//       if (convex_a->intersect(convex_b)) {
//         return true;
//       }
//     }
//   }

//   return false;
// }

/**
 * Creates a representation of all curves as linked lists between evaluated points.
 * Finds all intersections between curves and adds them as new points to the curves. These new
 * points contain references to the matching new point on the other curve. We're essentially
 * turning the splines into linked lists of points into linked meshes of points connected by common
 * pairs of intersection points.
 */
static std::list<CurveHull *> splineIntersectAll(std::vector<std::vector<Spline *>> splines)
{
  SCOPED_TIMER("Spline Intersect All");
  // TODO : Create Binary search tree from bounds (Needs some kind of clumping heuristics.)
  // TODO : Option to allow self intersection in exchange for brute force collision algorithm

  std::list<CurveHull *> results;
  std::list<CurveHull *> frontier;

  // Calculate bounding box for each individual curve.
  std::vector<Rect> bounds;
  // std::vector<std::vector<ConvexHull *>> convex_hulls;

  typedef boost::geometry::model::box<
      boost::geometry::model::point<float, 2UL, boost::geometry::cs::cartesian>>
      Box;
  typedef std::pair<Box, CurveHull *> BoxHullPair;
  // create the R* variant of the rtree (Shouldn't plain R version be better due to lower insert
  // cost?)
  bgi::rtree<BoxHullPair, bgi::rstar<8>> rtree;

  int c = splines.size();
  for (int i = 0; i < c; i++) {
    for (Spline *_spline : splines[i]) {
      float3 min = float3(INFINITY, INFINITY, INFINITY);
      float3 max = float3(-INFINITY, -INFINITY, -INFINITY);
      _spline->bounds_min_max(min, max, true);
      Rect bound = {(float2)min, (float2)max};
      bounds.push_back(bound);

      frontier.push_back(points_to_hull(_spline, i));
      CurveHull *last = frontier.back();
      int counter = count_containing_curves(last->curve_id == 0 ? splines[1] : splines[0],
                                            last->position);
      last->is_inside = counter % 2 == 1;

      //convex_hulls.push_back(generate_convex_hulls(last));

      // insert new value
      rtree.insert(std::make_pair(bound.Box(), last));
    }
  }

  // print("Full size " + str(rtree.size()) + " vs " + str(bounds.size()));

  // Boost takes longer, but feels faster?! Overall doesn't seem worth it, unless I can reduce the
  // actual line segments
  // TODO : Curve Hull binary search
  // Iterate CurveHull, parse into segments of multiple of resolution, but at least 4. Then build
  // binary tree using adjacent rects

  c = frontier.size();
  int pass = 0;
  int i = -1;  // refers to bounds array
  // loop through all possible curve segments.
  while (frontier.size() > 0) {
    i++;
    CurveHull *primary_curve = frontier.front();
    frontier.pop_front();

    // find values intersecting some area defined by a box
    std::vector<BoxHullPair> result_s;
    rtree.query(bgi::intersects(bounds[i].Box()), std::back_inserter(result_s));

    // print(convex_hulls.size());

    // int j = i;
    // for (BoxHullPair v : result_s) {
    //   j++;
    //   CurveHull *curve_hull_b = v.second;

    //   pass++;

    //   if (primary_curve->curve_id <= curve_hull_b->curve_id) {
    //     //||
    //     //!(bounds[i].intersect(Rect(float2(v.first.min_corner().get<0>(),v.first.min_corner().get<1>()),float2(v.first.max_corner().get<0>(),v.first.max_corner().get<1>()))))
    //     continue;
    //   }

    //   CurveHull *current_curve_hull_a = primary_curve;
    //   CurveHull *next_curve_hull_a = primary_curve->next;
    //   // loop through a's curve segments
    //   while (current_curve_hull_a != nullptr && next_curve_hull_a != nullptr &&
    //          current_curve_hull_a->pass != pass) {
    //     current_curve_hull_a->pass = pass;

    //     CurveHull *current_curve_hull_b = curve_hull_b;
    //     CurveHull *next_curve_hull_b = curve_hull_b->next;
    //     // loop through b's curve segments
    //     do {
    //       float2 intersection = linesIntersect(current_curve_hull_a->position,
    //                                            next_curve_hull_a->position,
    //                                            current_curve_hull_b->position,
    //                                            next_curve_hull_b->position);

    //       if (intersection[0] > 0) {  // Found intersection
    //         auto intersection_a = new
    //         Intersection{v : intersection[0], hull : current_curve_hull_a};
    //         auto intersection_b = new Intersection{
    //           v : intersection[1],
    //           target : intersection_a,
    //           hull : current_curve_hull_b
    //         };
    //         intersection_a->target = intersection_b;
    //         current_curve_hull_a->intersections.push_back(intersection_a);
    //         current_curve_hull_b->intersections.push_back(intersection_b);
    //       }

    //       current_curve_hull_b = current_curve_hull_b->next;
    //       next_curve_hull_b = current_curve_hull_b->next;
    //     } while (next_curve_hull_b != nullptr && current_curve_hull_b->curve_i != 0);

    //     current_curve_hull_a = next_curve_hull_a;
    //     next_curve_hull_a = current_curve_hull_a->next;
    //   }
    // }

    int j = i;  // refers to bounds array
    for (CurveHull *curve_hull_b : frontier) {
      j++;
      pass++;

      if (primary_curve->curve_id == curve_hull_b->curve_id)
        continue;

      // if (!any_collide(convex_hulls[i], convex_hulls[j])) {
      //   print("no collision");
      //   continue;
      // }
      if (primary_curve->curve_id == curve_hull_b->curve_id ||
      !(bounds[i].intersect(bounds[j]))) {
        continue;
      }

      CurveHull *current_curve_hull_a = primary_curve;
      CurveHull *next_curve_hull_a = primary_curve->next;
      // loop through a's curve segments
      while (current_curve_hull_a != nullptr && next_curve_hull_a != nullptr &&
             current_curve_hull_a->pass != pass) {
        current_curve_hull_a->pass = pass;

        CurveHull *current_curve_hull_b = curve_hull_b;
        CurveHull *next_curve_hull_b = curve_hull_b->next;
        // loop through b's curve segments
        do {
          float2 intersection = linesIntersect((float2)current_curve_hull_a->position,
                                               (float2)next_curve_hull_a->position,
                                               (float2)current_curve_hull_b->position,
                                               (float2)next_curve_hull_b->position);

          if (intersection[0] > 0) {  // Found intersection
            auto intersection_a = new
            Intersection{v : intersection[0], hull : current_curve_hull_a};
            auto intersection_b = new Intersection{
              v : intersection[1],
              target : intersection_a,
              hull : current_curve_hull_b
            };
            intersection_a->target = intersection_b;
            current_curve_hull_a->intersections.push_back(intersection_a);
            current_curve_hull_b->intersections.push_back(intersection_b);
          }

          current_curve_hull_b = current_curve_hull_b->next;
          next_curve_hull_b = current_curve_hull_b->next;
        } while (next_curve_hull_b != nullptr && current_curve_hull_b->curve_i != 0);

        current_curve_hull_a = next_curve_hull_a;
        next_curve_hull_a = current_curve_hull_a->next;
      }
    }

    results.push_back(primary_curve);
  }

  // Now that all intersections are stored in their relevant control points, sort them by v and
  // create new "virtual" CurveHull points for each intersection.
  for (CurveHull *current_point : results) {
    CurveHull *start = current_point;
    float curve_length = current_point->spline->evaluated_points_size() +
                         (current_point->spline->is_cyclic() ? 0 : -1);
    // Then create all intersection points as actual points. We didn't do this earlier to avoid
    // extra intersection calculations with the new line segments.
    do {
      std::vector<Intersection *> &intersections = current_point->intersections;
      // Sort the intersections by progress along line (v).
      std::sort(intersections.begin(), intersections.end(), [](Intersection *a, Intersection *b) {
        return a->v < b->v;
      });

      CurveHull *original_current_point = current_point;
      CurveHull *original_next_point = current_point->next;
      for (Intersection *intersection : intersections) {
        CurveHull *new_point = new CurveHull{
            .curve_id = current_point->curve_id,
            .spline = current_point->spline,
            .curve_i = current_point->curve_i,
            .v = original_current_point->v + intersection->v / curve_length,
            .position = math::interpolate(original_current_point->position, original_next_point->position, intersection->v),
            .next = current_point->next,
            .prev = current_point,
            .intersection = intersection};
        intersection->hull = new_point;
        current_point->next->prev = new_point;
        current_point->next = new_point;
        current_point = new_point;
      }

      current_point = current_point->next;
    } while (current_point != nullptr && current_point != start);
  }

  return results;
}

/**
 * Given a bezier curve segment described by 2 points and 2 handles, generate an extra point at v
 * along the curve and return both the new point's position and handles, as well as the modified
 * handles of the initial points.
 * This is the same as #BezierSpline::calculate_bezier_segment_insertion, but it doesn't need the
 * points to exist on the curve yet.
 */
static BezierSpline::InsertResult calculate_bezier_segment_insertion(
    float3 pos_prev, float3 handle_prev, float3 pos_next, float3 handle_next, float v)
{
  const float3 center_point = math::interpolate(handle_prev, handle_next, v);

  BezierSpline::InsertResult result;
  result.handle_prev = math::interpolate(pos_prev, handle_prev, v);
  result.handle_next = math::interpolate(handle_next, pos_next, v);
  result.left_handle = math::interpolate(result.handle_prev, center_point, v);
  result.right_handle = math::interpolate(center_point, result.handle_next, v);
  result.position = math::interpolate(result.left_handle, result.right_handle, v);
  return result;
}

/**
 * Assign bezier point data to a particular index in a #BezierSpline.
 */
static void set_bezier_point(Spline *spline,
                             int index,
                             float3 position,
                             float3 handle_l,
                             float3 handle_r,
                             float radius,
                             float tilt)
{
  BezierSpline *bezierSpline = (BezierSpline *)spline;
  bezierSpline->positions()[index] = position;
  bezierSpline->handle_positions_left()[index] = handle_l;
  bezierSpline->handle_positions_right()[index] = handle_r;
  bezierSpline->radii()[index] = radius;
  bezierSpline->tilts()[index] = tilt;
  bezierSpline->handle_types_left()[index] = HandleType::Free;
  bezierSpline->handle_types_right()[index] = HandleType::Free;
}

/**
 * Assign spline point data to a particular index in a #Spline of any type.
 */
static void set_spline_point(Spline *spline, int index, float3 position, float radius, float tilt)
{
  BezierSpline *bezierSpline = (BezierSpline *)spline;
  bezierSpline->positions()[index] = position;
  bezierSpline->radii()[index] = radius;
  bezierSpline->tilts()[index] = tilt;
}

/**
 * Checks if the current point is a control point.
 * If so, inserts control point data at the specified index.
 */
static float process_control_point(CurveHull *current_point,
                                   Spline *result,
                                   int result_index,
                                   Spline::Type type,
                                   float last_intersection_offset = 0)
{
  int resolution = current_point->spline->evaluated_points_size() /
                   current_point->spline->segments_size();

  // control point index + v along segment
  float current_absolute_v = current_point->curve_i / (float)(resolution);
  int current_control_point = (int)(current_absolute_v);

  if (type == Spline::Type::Bezier) {
    if ((current_point->curve_i % resolution) != 0) {
      return last_intersection_offset;
    }
    BezierSpline *current_spline = (BezierSpline *)current_point->spline;
    BezierSpline *result_spline = (BezierSpline *)result;

    set_bezier_point(
        result,
        result_index,
        current_spline->evaluated_positions()[current_point->curve_i],
        math::interpolate(current_spline->handle_positions_left()[current_control_point],
                            current_spline->positions()[current_control_point],
                            last_intersection_offset),
        current_spline->handle_positions_right()[current_control_point],
        current_spline->radii()[current_control_point],
        current_spline->tilts()[current_control_point]);
  }
  else {
    Spline *current_spline = current_point->spline;

    int current_segment_end = (current_control_point + 1) % current_spline->size();
    float current_relative_v = fmod(current_absolute_v, 1);
    MutableSpan<float> radii = current_spline->radii();
    MutableSpan<float> tilts = current_spline->tilts();
    set_spline_point(
        result,
        result_index,
        current_spline->evaluated_positions()[current_point->curve_i],
        // We mix, because if this is a bezier curve, that's forcefully converted to poly, it won't
        // have evaluated radii
        mix2(current_relative_v, radii[current_control_point], radii[current_segment_end]),
        mix2(current_relative_v, tilts[current_control_point], tilts[current_segment_end]));
  }

  return 0;
}

/**
 * Process a #CurveHull point that is an intersection.
 * Inserts control point data at the specified index in the result spline.
 */
static void process_intersection(CurveHull *current_point,
                                 Spline *result,
                                 int result_index,
                                 float &last_intersection_offset,
                                 Spline::Type type)
{
  bool first = result_index == 0;  // is this the first point in the newly generated spline?
  Intersection *current_intersection = current_point->intersection;

  CurveHull *other_point = current_intersection->target->hull;

  if (type == Spline::Type::Bezier) {
    BezierSpline *bezier_result = (BezierSpline *)result;

    BezierSpline *current_spline = (BezierSpline *)current_point->spline;
    BezierSpline *other_spline = (BezierSpline *)other_point->spline;

    MutableSpan<float3> current_positions = current_spline->positions();
    int current_resolution = current_spline->resolution();
    float current_absolute_v = current_point->curve_i / (float)(current_resolution);
    int current_control_point = (int)(current_absolute_v);
    float current_relative_v = fmod(current_absolute_v, 1);
    int current_segment_end = (current_control_point + 1) % current_spline->size();
    float v_cur = (current_relative_v + current_intersection->v / current_resolution -
                   last_intersection_offset) /
                  (1 - last_intersection_offset);

    float3 last_pos = first ? current_positions[current_control_point] :
                              result->positions()[result_index - 1];
    float3 last_handle = first ? current_spline->handle_positions_right()[current_control_point] :
                                 bezier_result->handle_positions_right()[result_index - 1];
    BezierSpline::InsertResult insert_current = calculate_bezier_segment_insertion(
        last_pos,
        last_handle,
        current_positions[current_segment_end],
        math::interpolate(current_spline->handle_positions_left()[current_segment_end],
                            current_positions[current_segment_end],
                            last_intersection_offset),
        v_cur);

    if (!first) {
      bezier_result->handle_positions_right()[result_index - 1] = insert_current.handle_prev;
    }

    MutableSpan<float3> other_positions = other_spline->positions();
    int next_resolution = other_spline->resolution();
    float next_absolute_v = other_point->curve_i / (float)(next_resolution);
    int next_control_point = (int)(next_absolute_v);
    float next_relative_v = fmod(next_absolute_v, 1);
    int next_segment_end = (next_control_point + 1) % other_spline->size();
    last_intersection_offset = next_relative_v + current_intersection->target->v / next_resolution;

    BezierSpline::InsertResult insert_next = calculate_bezier_segment_insertion(
        other_positions[next_control_point],
        other_spline->handle_positions_right()[next_control_point],
        other_positions[next_segment_end],
        other_spline->handle_positions_left()[next_segment_end],
        last_intersection_offset);

    float3 real_pos =
        current_intersection->hull->position;  // the intersection point on the EVALUATED geometry

    MutableSpan<float> radii = current_spline->radii();
    MutableSpan<float> tilts = current_spline->tilts();
    set_bezier_point(
        result,
        result_index,
        real_pos,
        insert_current.left_handle,  //+ real_pos - bezier_pos_current
        insert_next.right_handle,    // + real_pos - bezier_pos_next
        mix2(next_relative_v, radii[current_control_point], radii[current_segment_end]),
        mix2(next_relative_v, tilts[current_control_point], tilts[current_segment_end]));
  }
  else {
    Spline *current_spline = current_point->spline;

    int current_resolution = current_spline->evaluated_points_size() /
                             current_point->spline->segments_size();
    float current_absolute_v = current_point->curve_i / (float)(current_resolution);
    int current_control_point = (int)(current_absolute_v);
    int current_segment_end = (current_control_point + 1) % current_spline->size();

    MutableSpan<float> radii = current_spline->radii();
    MutableSpan<float> tilts = current_spline->tilts();
    float3 real_pos = current_intersection->hull->position;
    set_spline_point(
        result,
        result_index,
        real_pos,
        mix2(current_intersection->v, radii[current_control_point], radii[current_segment_end]),
        mix2(current_intersection->v, tilts[current_control_point], tilts[current_segment_end]));
  }
}

/**
 * Given a certain starting point, move along the #CurveHull and store each point in the result
 * list. Turns at each intersection. Returns an empty list if the trace failed to find a valid
 * path.
 */
static std::list<CurveHull *> do_trace(CurveHull *start_point,
                                       std::list<CurveHull *> &frontier,
                                       Spline::Type return_type)
{
  std::list<CurveHull *> result;

  CurveHull *current_point = start_point;
  CurveHull *next_point = current_point->next;

  int resolution = return_type == Spline::Type::Bezier ?
                       ((BezierSpline *)start_point->spline)->resolution() :
                       1;

  if (next_point != nullptr) {
    while (current_point != nullptr && (current_point->pass != -2)) {
      if (current_point->intersection == nullptr && (current_point->curve_i % resolution) != 0) {
        current_point = current_point->next;
        continue;  // only keep control and intersection points.
      }
      result.push_back(current_point);

      current_point->pass = -2;

      if (current_point->intersection != nullptr)  // current point is a virtual intersection point
      {
        CurveHull *next_possible_frontier = current_point->next_intersection();
        if (next_possible_frontier != nullptr) {
          frontier.push_back(next_possible_frontier->intersection->target->hull);
        }

        current_point = current_point->intersection->target->hull;
        current_point->pass = -2;
        resolution = return_type == Spline::Type::Bezier ?
                         ((BezierSpline *)current_point->spline)->resolution() :
                         1;
      }
      current_point = current_point->next;
    }

    if (current_point != nullptr) {
      if ((current_point != start_point) &&
          (current_point->intersection == nullptr ||
           current_point->intersection->target->hull != start_point)) {
        // skip, didn't reach a valid target. Such as the part between two processed intersections.
        return {};
      }
    }
  }

  return result;
}

/**
 * After knowing a trace is viable, create a #Spline of type `return_type` from it.
 */
static std::unique_ptr<Spline> trace_to_spline(std::list<CurveHull *> trace,
                                               Spline::Type return_type)
{
  std::unique_ptr<Spline> result = return_type == Spline::Type::Bezier ?
                                       (std::unique_ptr<Spline>)std::make_unique<BezierSpline>() :
                                       (std::unique_ptr<Spline>)std::make_unique<PolySpline>();
  result->resize(trace.size());

  float last_intersection_offset = 0;
  int i = 0;
  for (CurveHull *current_point : trace) {
    if (current_point->intersection != nullptr)  // current point is a virtual intersection point
    {
      process_intersection(current_point, result.get(), i, last_intersection_offset, return_type);
    }
    else {
      // This function skips non control points by itself.
      last_intersection_offset = process_control_point(
          current_point, result.get(), i, return_type, last_intersection_offset);
    }
    i++;
  }

  CurveHull *start_point = trace.front();

  if (return_type == Spline::Type::Bezier)  // if cyclical result
  {
    BezierSpline *bezier_result = (BezierSpline *)result.get();
    // if cyclic, make sure the initial point's left handle scales with any previous intersection.
    // (If not cyclic, those handles aren't used)
    if (start_point->intersection != nullptr) {
      int segment_count = start_point->spline->segments_size();

      float current_relative_v = fmod((start_point->v * segment_count), 1.0f);

      // both handles are already resized to the following sizes
      float expected_distance_right = 1 - last_intersection_offset;
      float expected_distance_left = current_relative_v;

      // but because of this new intersection, they will need to be this size instead
      float proportional_distance_right = (1 - current_relative_v) / expected_distance_right;
      float proportional_distance_left = last_intersection_offset / expected_distance_left;

      int back_i = bezier_result->size() - 1;
      bezier_result->handle_positions_left()[0] = math::interpolate(
          bezier_result->handle_positions_left()[0],
          bezier_result->positions()[0],
          proportional_distance_left);
      bezier_result->handle_positions_right()[back_i] = math::interpolate(
          bezier_result->handle_positions_right()[back_i],
          bezier_result->positions()[back_i],
          proportional_distance_right);
    }
    else {
      bezier_result->handle_positions_left()[0] = math::interpolate(
          bezier_result->handle_positions_left()[0],
          bezier_result->positions()[0],
          last_intersection_offset);
    }
  }

  return result;
}

/**
 * Walk along the #CurveHull and generate new #Splines from each valid path.
 */
static void trace_hull(std::vector<std::unique_ptr<Spline>> &results,
                       CurveHull *path,
                       BoolOperationType type,
                       Spline::Type return_type)
{
  std::list<CurveHull *> frontier;
  frontier.push_back(path);

  while (frontier.size() > 0) {
    CurveHull *start_point = frontier.front();
    frontier.pop_front();

    bool is_inside = (type == BoolOperationType::AND) ? start_point->is_inside == 0 :
                                                        start_point->is_inside == 1;

    if (is_inside) {
      start_point = start_point->next_intersection();
      if (start_point == nullptr) {
        continue;
      }
      start_point = start_point->intersection->target->hull;
    }

    if (start_point->pass == -2) {
      continue;
    }

    std::list<CurveHull *> trace = do_trace(start_point, frontier, return_type);

    if (trace.size() > 0) {
      results.push_back(trace_to_spline(trace, return_type));
    }
  }
}

/**
 * Do the evaluated points form a clockwise curve?
 * #Spline must have at least one point.
 */
static bool is_curve_clockwise(Spline *spline)
{
  Span<float3> points = spline->evaluated_positions();

  float dotsum = 0;  // calculates 2x enclosed area. If negative, curve is counter clockwise
  for (int i = 0; i < points.size() - 1; i++) {
    dotsum += (points[i + 1].x - points[i].x) * (points[i + 1].y + points[i].y);
  }

  const float3 start = points.first();
  const float3 end = points.last();
  dotsum += (start.x - end.x) * (start.y + end.y);

  return dotsum >= 0;
}

/**
 * Extract #Spline Pointers and put them into a resizable vector.
 * Toss all splines with length 0 for one less border case.
 */
static std::vector<Spline *> _unrwap_splines(blender::Span<SplinePtr> s)
{
  std::vector<Spline *> result = {};
  for (const SplinePtr &spline : s) {
    if (spline.get()->size() > 0)
      result.push_back(spline.get());
  }
  return result;
}

/**
 * If a curve is a bezier curve, switch all control points to mode "free".
 */
static void free_handles(Spline *spline)
{
  if (spline->type() == Spline::Type::Bezier) {
    BezierSpline *bezier_spline = ((BezierSpline *)spline);
    MutableSpan<HandleType> handles_l_t = bezier_spline->handle_types_left();
    MutableSpan<HandleType> handles_r_t = bezier_spline->handle_types_right();
    bezier_spline->ensure_auto_handles();

    int size = handles_l_t.size();
    for (int i = 0; i < size; i++) {
      handles_l_t[i] = HandleType::Free;
      handles_r_t[i] = HandleType::Free;
    }

    spline->mark_cache_invalid();
  }
}

/**
 * Generate one or more new curves from 2 existing sets of curves.
 * These curves must not self intersect.
 * The general idea is to follow one of the curves and copy the control points until an
 * intersection is found. At which point we know the other curve is bigger, so we switch to the
 * other one. We do this for each intersection until we reach our initial position.
 */
static std::unique_ptr<CurveEval> generate_boolean_shapes(const CurveEval *a,
                                                          const CurveEval *b,
                                                          BoolOperationType type,
                                                          int resolution)
{
  std::unique_ptr<CurveEval> result = std::make_unique<CurveEval>();

  // SplinePtr is a unique ptr. We're not modifying those, we just want their
  // data. So we cast to Spline* from the get-go.
  std::vector<Spline *> splines_a = _unrwap_splines(a->splines());
  std::vector<Spline *> splines_b = _unrwap_splines(b->splines());
  std::vector<std::vector<Spline *>> splines = {splines_a, splines_b};

  Spline::Type result_type = Spline::Type::Bezier;

  for (Span<Spline *> _splines : splines) {
    for (Spline *spline : _splines) {
      if (spline->type() != Spline::Type::Bezier) {
        result_type = Spline::Type::Poly;
      }
    }
  }

  /// First, we need to prepare the curves. Outside curves must go clockwise, inside
  /// counterclockwise. Operators may flip them yet again.

  for (Span<Spline *> _splines : splines) {
    for (Spline *spline : _splines) {
      // Auto Handles are weird. Make sure all handles are Free.
      free_handles(spline);

      // All curves must be clockwise, unless the algorithm turns them explicitly.
      bool is_clockwise = is_curve_clockwise(spline);

      // Outside curves must go clockwise.
      // Inside curves must go counterclockwise.
      int counter = count_containing_curves(_splines, spline->positions()[0], spline);
      if (counter % 2 == (is_clockwise ? 1 : 0)) {
        // curve represents negative space, must go counterclockwise.
        spline->reverse();
      }
    }
  }

  // Some operations need one curve to go clockwise, and the other counter-clockwise.
  if (type == BoolOperationType::SUB) {
    for (Spline *spline : splines_b) {
      spline->reverse();
    }
  }

  /// Then, find all intersections
  std::list<CurveHull *> paths = splineIntersectAll(splines);

  int length = 0;

  for (CurveHull *path : paths) {
    std::vector<std::unique_ptr<Spline>> results;
    if (path->next_intersection() == nullptr)  // Path does not intersect with anything.
    {
      // count the number of curves that contain this one. Depending on the mode, and if it's
      // contained in the other curve, toss it, or keep it as-is.
      // Eg if the mode is Union, and curve B is fully inside curve A, it won't do anything,
      // so we discard it.
      int counter = count_containing_curves(path->curve_id == 0 ? splines_b : splines_a,
                                            path->position);
      if (type == BoolOperationType::SUB) {
        if (counter % 2 == (path->curve_id == 0 ? 0 : 1)) {
          results.push_back(path->spline->copy());
        }
      }
      else if (type == BoolOperationType::OR) {
        if (counter % 2 == 0) {
          results.push_back(path->spline->copy());
        }
      }
      else if (type == BoolOperationType::AND) {
        if (counter % 2 == 1) {
          results.push_back(path->spline->copy());
        }
      }
    }
    else {
      trace_hull(results, path, type, result_type);
    }

    /// Now turn our individual spline pointers into their final unique ptr form.

    for (std::unique_ptr<Spline> &trace : results) {
      if (trace->type() == Spline::Type::Bezier) {
        BezierSpline *result_spline = (BezierSpline *)trace.get();
        // BUG : Why does removing this cause SegFaults?
        result_spline->set_resolution(resolution);
      }
      trace.get()->set_cyclic(true);
      result->add_spline(std::move(trace));
      length += result->splines().size();
    }
  }

  // Remove all CurveHull and Intersection helper structs. These were created with new().
  for (CurveHull *path : paths) {
    path->remove();
  }

  result->attributes.reallocate(length);

  return result;
}

/**
 * This is called whenever the node needs to update, such as when inputs change.
 */
static void geo_node_curve_bool_exec(GeoNodeExecParams params)
{
  SCOPED_TIMER("Curve Boolean");

  int resolution = params.extract_input<int>("Int. Resolution");
  GeometrySet geometry_set_a = params.extract_input<GeometrySet>("Base Curve");
  //CurveComponent& curve_set_a = geometry_set_a.get_component_for_write<CurveComponent>();//bke::geometry_set_realize_instances(curve_set_a);
  std::vector<const CurveEval *> curves_b;

  Vector<GeometrySet> geometry_sets = params.extract_multi_input<GeometrySet>("Curves");

  if (!(geometry_set_a.has_curve())) {
    print("Nope");
    params.set_output("Curve", geometry_set_a);
    return;
  }

  const bNode &node = params.node();
  const BoolOperationType data_type = static_cast<BoolOperationType>(node.custom2);

  int real_count = 0;
  // Repeat the bool operation for each element in the operand `Curves` input.
  const CurveEval *primary_curve = geometry_set_a.get_curve_for_read();
  std::unique_ptr<CurveEval> result_curve;
  for (const GeometrySet &set_group : geometry_sets) {
    // const CurveComponent* curve_in = set_group.get_component_for_read<CurveComponent>();
    if (set_group.has_curve()) {
      result_curve = generate_boolean_shapes(
          primary_curve, set_group.get_curve_for_read(), data_type, resolution);

      primary_curve = result_curve.get();

      real_count++;
    }
  }

  if (real_count == 0) {
    params.set_output("Curve", geometry_set_a);
  }
  else {
     params.set_output("Curve", GeometrySet::create_with_curve(result_curve.release()));
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
