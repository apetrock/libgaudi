// Explicit template instantiations for gaudi::arp library
// This file compiles the commonly used template instantiations once
// to avoid duplicate compilation in every object file

#include "gaudi/arp/morton.hpp"
#include "gaudi/arp/hash_tree.hpp"

namespace gaudi {
namespace arp {

// ============================================================================
// Legacy explicit instantiations (flat Vec3View)
// ============================================================================

// Explicit instantiations for calc_com<N> (from morton.hpp)
template std::vector<MassPoint> calc_com<1>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<2>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<3>(const std::vector<vec3>&);

// Explicit instantiations for calc_extents<N> (from morton.hpp)
template std::vector<ext::extents_t> calc_extents<1>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<2>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<3>(const std::vector<vec3>&);

// Explicit instantiations for make_points<N> (from hash_tree.hpp)
template TreeResult<vec3> make_points<1>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<2>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<3>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// Explicit instantiations for make_bvh<N> (from hash_tree.hpp)
template TreeResult<ext::extents_t> make_bvh<1>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<2>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<3>(const std::vector<vec3>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// ============================================================================
// New explicit instantiations for SimplexView (type-based templating)
// ============================================================================

// Type aliases for common SimplexView types
using edge_simplex_view = simplex_view<2, std::vector<vec3>>;
using tri_simplex_view = simplex_view<3, std::vector<vec3>>;
using point_simplex_view = simplex_view<1, std::vector<vec3>>;

using edge_psv = permuted_simplex_view<2, std::vector<vec3>, std::vector<index_t>>;
using tri_psv = permuted_simplex_view<3, std::vector<vec3>, std::vector<index_t>>;
using point_psv = permuted_simplex_view<1, std::vector<vec3>, std::vector<index_t>>;

// calc_com for SimplexView types
template std::vector<MassPoint> calc_com<edge_simplex_view>(const edge_simplex_view&);
template std::vector<MassPoint> calc_com<tri_simplex_view>(const tri_simplex_view&);
template std::vector<MassPoint> calc_com<point_simplex_view>(const point_simplex_view&);
template std::vector<MassPoint> calc_com<edge_psv>(const edge_psv&);
template std::vector<MassPoint> calc_com<tri_psv>(const tri_psv&);
template std::vector<MassPoint> calc_com<point_psv>(const point_psv&);

// calc_extents for SimplexView types
template std::vector<ext::extents_t> calc_extents<edge_simplex_view>(const edge_simplex_view&);
template std::vector<ext::extents_t> calc_extents<tri_simplex_view>(const tri_simplex_view&);
template std::vector<ext::extents_t> calc_extents<point_simplex_view>(const point_simplex_view&);
template std::vector<ext::extents_t> calc_extents<edge_psv>(const edge_psv&);
template std::vector<ext::extents_t> calc_extents<tri_psv>(const tri_psv&);
template std::vector<ext::extents_t> calc_extents<point_psv>(const point_psv&);

// make_hash for SimplexView types
template std::tuple<std::vector<morton_t>, std::vector<index_t>,
                    std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash<morton_t, edge_psv>(const edge_psv&);
template std::tuple<std::vector<morton_t>, std::vector<index_t>,
                    std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash<morton_t, tri_psv>(const tri_psv&);
template std::tuple<std::vector<morton_t>, std::vector<index_t>,
                    std::vector<radix_tree_node>, std::vector<radix_tree_node>>
make_hash<morton_t, point_psv>(const point_psv&);

// make_bvh for SimplexView types
template TreeResult<ext::extents_t> make_bvh<edge_psv>(const edge_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<tri_psv>(const tri_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<point_psv>(const point_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// make_points for SimplexView types
template TreeResult<vec3> make_points<edge_psv>(const edge_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<tri_psv>(const tri_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<point_psv>(const point_psv&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// bvh_tree_t for SimplexView types
template class bvh_tree_t<edge_psv>;
template class bvh_tree_t<tri_psv>;
template class bvh_tree_t<point_psv>;

} // namespace arp
} // namespace gaudi 