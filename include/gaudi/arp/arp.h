#ifndef __ARP__
#define __ARP__
/*
https://en.wikipedia.org/wiki/Jean_Arp
Hans Peter Wilhelm Arp (16 September 1886 – 7 June 1966), better known as Jean
Arp in English, was a German-French sculptor, painter, and poet. He was known as
a Dadaist and an abstract artist.
*/

// ---------------------------------------------------------------------------
// Backend selection (module export surface)
// ---------------------------------------------------------------------------
// USE_HASH selects the acceleration structure backing the arp proximity
// queries:
//   defined   -> 64-bit Morton / radix BVH (hash_tree.hpp)
//   undefined -> legacy half-space AABB tree (aabb.hpp)       [current]
//
// Flip it here, in one place, to validate the Morton path against the legacy
// tree (e.g. when a 64-bit hash collision is suspected). The fast-summation /
// Barnes-Hut path always needs the radix tree, so hash_tree.hpp is included
// unconditionally; the legacy header is pulled in as the fallback backend.

#define USE_HASH

#include "gaudi/arp/hash_tree.hpp"

#if !defined(USE_HASH)
#include "gaudi/arp/aabb.hpp"
#endif

// ---------------------------------------------------------------------------
// Backend-neutral proximity / FMM tree aliases (the export surface).
// ---------------------------------------------------------------------------
// Every consumer talks to arp::T1/T2/T3; USE_HASH chooses which acceleration
// structure backs them. Both backends expose the same canonical representation
// (arp::tree_view) plus create/update/find_nearest/find_neighbors/get_nearest
// and the leaf_simplex/get_com/bvh_ accessors the calder FMM path relies on, so
// flipping the toggle requires no consumer changes.
namespace gaudi {
namespace arp {
#if defined(USE_HASH)
using T1 = bvh_tree<1, morton_t>;
using T2 = bvh_tree<2, morton_t>;
using T3 = bvh_tree<3, morton_t>;
#else
using T1 = aabb_tree<1>;
using T2 = aabb_tree<2>;
using T3 = aabb_tree<3>;
#endif
} // namespace arp
} // namespace gaudi

#endif
