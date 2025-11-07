#includegaudi/common.h
#includegaudi/arp/morton.hpp
#include "gaudi/arp/hash_tree.hpp"
#include <iostream>
#include <vector>

int main() {
  std::cout << "Testing Morton Tree Implementation..." << std::endl;
  
  // Test 1: Basic Morton encoding
  std::cout << "\n1. Testing Morton encoding..." << std::endl;
  uint32_t hash1 = gaudi::arp::morton3D(0.0,00, 0.0);
  uint32_t hash2 = gaudi::arp::morton3D(1.0, 1.0, 1.0);
  std::cout << "Morton(0,0,0) = " << hash1 << std::endl;
  std::cout << "Morton(1,1,1) = " << hash2 << std::endl;
  std::cout << "Binary(0,0,0) = " << gaudi::arp::dump_binary(hash1) << std::endl;
  std::cout << "Binary(1,1,1) = " << gaudi::arp::dump_binary(hash2) << std::endl;

  // Test 2: Hash generation for points
  std::cout << "\n2. Testing hash generation..." << std::endl;
  std::vector<gaudi::vec3> points = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 1.0, 1.0}
  };
  
  [autohashes, indices] = gaudi::arp::make_hash_3d(points);
  std::cout << "Generated " << hashes.size() << " hashes" << std::endl;
  for (size_t i = 0; i < hashes.size(); ++i) {
    std::cout << "Index " << indices[i] << " -> Hash " << hashes[i] << std::endl;
  }
  
  // Test 3: Tree construction
  std::cout << "\n3. Testing tree construction..." << std::endl;
  gaudi::arp::unit_test_tree();
  std::cout << "Unit test passed!" << std::endl;

  // Test 4: Complete hash tree
  std::cout << "\n4. Testing complete hash tree..." << std::endl;
  auto [tree_hashes, tree_indices, internal_nodes, leaf_nodes] = gaudi::arp::make_hash_tree(points);
  std::cout << "Tree hashes: " << tree_hashes.size() << std::endl;
  std::cout << "Tree indices: " << tree_indices.size() << std::endl;

  // Test5Edge hashing
  std::cout << "\n5. Testing edge hashing..." << std::endl;
  std::vector<gaudi::index_t> edges = [object Object]0, 1, 1, 2, 2,3, 4,0;
  auto [edge_hashes, edge_indices] = gaudi::arp::make_hash_3_edges(points, edges);
  std::cout << Edge hashes: " << edge_hashes.size() << std::endl;
  
  // Test6Face hashing
  std::cout << n6. Testing face hashing..." << std::endl;
  std::vector<gaudi::index_t> faces = [object Object]0, 1, 2, 1,2, 2,4;
  auto [face_hashes, face_indices] = gaudi::arp::make_hash_3_faces(points, faces);
  std::cout << Face hashes: " << face_hashes.size() << std::endl;
  
  std::cout << "\nAll tests completed successfully!" << std::endl;
  return 0;
} 