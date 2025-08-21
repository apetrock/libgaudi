/**
 * TypeScript interface for morton_tree WASM module
 * Generated automatically - do not edit manually
 */

export interface MortonTreeTestClass {
  test_morton_encoding(x: number, y: number, z: number): number;
  get_binary_string(hash: number): string;
  get_hash_count(): number;
  get_hash(index: number): number;
  get_index(index: number): number;
  test_tree_construction(): boolean;
  test_complete_hash_tree(): boolean;
  test_edge_hashing(): boolean;
  test_face_hashing(): boolean;
  add_point(x: number, y: number, z: number): void;
  clear_points(): void;
  get_point_count(): number;
  generate_random_points(count: number): void;
  generate_grid_points(grid_size: number): void;
  generate_trefoil_knot(num_segments: number): void;
  generate_knot(k: number, num_segments: number): void;
  generate_random_knot(num_control_points: number, segments_per_chord: number): void;
  mk_hash_tree(): boolean;
  log_hierarchy(): void;
  log_bvh(): void;
  log_nearest(t: number): void;
  get_point_x(index: number): number;
  get_point_y(index: number): number;
  get_point_z(index: number): number;
  get_leaf_point_x(index: number): number;
  get_leaf_point_y(index: number): number;
  get_leaf_point_z(index: number): number;
  get_leaf_count(): number;
  
  // Visualization controls
  set_auto_visualize(enable: boolean): void;
  get_auto_visualize(): boolean;
  clear_visualizations(): void;
  draw_points_only(): void;
  draw_sorted_lines_only(): void;
  draw_all_visualizations(): void;
  
  // Logger API compatibility
  get_line_count(): number;
  get_point_count_logger(): number;
  get_line(index: number, start: Float64Array, end: Float64Array, color: Float64Array): void;
  get_point_logger(index: number, position: Float64Array, color: Float64Array): void;
}

export interface MortonTreeTestModule {
  MortonTreeTest: new () => MortonTreeTestClass;
}

/**
 * Morton tree test WASM module
 */
export const morton_treeEndpoints = {
  moduleName: 'morton_tree_test',
  functions: [
    { name: 'test_morton_encoding', returnType: 'number', parameters: [{ name: 'x', type: 'number' }, { name: 'y', type: 'number' }, { name: 'z', type: 'number' }] },
    { name: 'get_binary_string', returnType: 'string', parameters: [{ name: 'hash', type: 'number' }] },
    { name: 'get_hash_count', returnType: 'number', parameters: [] },
    { name: 'get_hash', returnType: 'number', parameters: [{ name: 'index', type: 'number' }] },
    { name: 'get_index', returnType: 'number', parameters: [{ name: 'index', type: 'number' }] },
    { name: 'test_tree_construction', returnType: 'boolean', parameters: [] },
    { name: 'test_complete_hash_tree', returnType: 'boolean', parameters: [] },
    { name: 'test_edge_hashing', returnType: 'boolean', parameters: [] },
    { name: 'test_face_hashing', returnType: 'boolean', parameters: [] },
    { name: 'add_point', returnType: 'void', parameters: [{ name: 'x', type: 'number' }, { name: 'y', type: 'number' }, { name: 'z', type: 'number' }] },
    { name: 'clear_points', returnType: 'void', parameters: [] },
    { name: 'get_point_count', returnType: 'number', parameters: [] },
    { name: 'generate_trefoil_knot', returnType: 'void', parameters: [{ name: 'num_segments', type: 'number' }] },
    { name: 'generate_knot', returnType: 'void', parameters: [{ name: 'k', type: 'number' }, { name: 'num_segments', type: 'number' }] },
    { name: 'generate_random_knot', returnType: 'void', parameters: [{ name: 'num_control_points', type: 'number' }, { name: 'segments_per_chord', type: 'number' }] }
  ]
} as const; 