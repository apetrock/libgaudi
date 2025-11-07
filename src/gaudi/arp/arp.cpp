// ARP Library Implementation
// This file contains all explicit template instantiations and library functions
// for the gaudi::arp namespace

#include "gaudi/arp/arp.h"
#include "gaudi/console_logger.hpp"

namespace gaudi {
namespace arp {

// Library version
const char* get_arp_library_version() {
    return "1.0.0";
}

// Library initialization
void initialize_arp_library() {
    console_logger::info << "Initializing ARP library v" << get_arp_library_version() << std::endl;
}

// Library cleanup
void cleanup_arp_library() {
    console_logger::info << "Cleaning up ARP library" << std::endl;
}

} // namespace arp
} // namespace gaudi

// Explicit template instantiations for the ARP library
// These are compiled once and linked by users who want pre-compiled versions

// From morton.hpp
namespace gaudi {
namespace arp {

// Explicit instantiations for calc_com<N>
template std::vector<MassPoint> calc_com<1>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<2>(const std::vector<vec3>&);
template std::vector<MassPoint> calc_com<3>(const std::vector<vec3>&);

// Explicit instantiations for calc_extents<N>
template std::vector<ext::extents_t> calc_extents<1>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<2>(const std::vector<vec3>&);
template std::vector<ext::extents_t> calc_extents<3>(const std::vector<vec3>&);

// From hash_tree.hpp
// Explicit instantiations for make_points<N>
template TreeResult<vec3> make_points<1>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<2>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<vec3> make_points<3>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// Explicit instantiations for make_bvh<N>
template TreeResult<ext::extents_t> make_bvh<1>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<2>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);
template TreeResult<ext::extents_t> make_bvh<3>(const std::vector<vec3>&, const std::vector<index_t>&, const std::vector<uint32_t>&, const std::vector<radix_tree_node>&, const std::vector<radix_tree_node>&);

// Explicit instantiations for make_hash_N<N>
template std::tuple<std::vector<uint32_t>, std::vector<index_t>, std::vector<radix_tree_node>, std::vector<radix_tree_node>> make_hash(const std::vector<vec3>&);
} // namespace arp
} // namespace gaudi 