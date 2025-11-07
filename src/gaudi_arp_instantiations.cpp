// Explicit template instantiations for gaudi::arp library
// This file compiles the commonly used template instantiations once
// to avoid duplicate compilation in every object file

#include "gaudi/arp/morton.hpp"
#include "gaudi/arp/hash_tree.hpp"

namespace gaudi {
namespace arp {

// Explicit instantiations for calc_com<N> (from morton.hpp)
template std::vector<MassPoint> calc_com<1>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<2>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<3>(const std::vector<vec3>&);

// Explicit instantiations for calc_extents<N> (from morton.hpp)
template std::vector<ext::extents_t> calc_extents<1>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<2>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<3>(const std::vector<vec3>&);

// Explicit instantiations for make_points<N> (from hash_tree.hpp)
template TreeResult<vec3> make_points<1>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<2>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<3>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// Explicit instantiations for make_bvh<N> (from hash_tree.hpp)
template TreeResult<ext::extents_t> make_bvh<1>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<2>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<3>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

} // namespace arp
} // namespace gaudi 