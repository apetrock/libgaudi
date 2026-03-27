export interface ModuleConfig {
  name: string;
  description: string;
  sources: string[];
  useLogger: boolean;
  exportName?: string;
  linkFlags?: string[];
  outputFiles: string[];
}

export const MODULE_CONFIGS: Record<string, ModuleConfig> = {
  'morton_tree': {
    name: 'morton_tree',
    description: 'Morton Tree Test Module',
    sources: ['src/morton_tree_test.cpp'],
    useLogger: true,
    exportName: 'MortonTreeTestModule',
    outputFiles: ['morton_tree_test.js']
  },
  /*
  'logger': {
    name: 'logger',
    description: 'Logger API',
    sources: ['src/gaudi_logger_test.cpp'],
    useLogger: true,
    exportName: 'GaudiLoggerTestModule',
    outputFiles: ['gaudi_logger_test.js']
  },
  */

  'rod_constraints': {
    name: 'rod_constraints',
    description: 'Rod Constraints Simulation',
    sources: ['src/rod_constraints_test.cpp'],
    useLogger: true,
    exportName: 'RodConstraintsTestModule',
    outputFiles: ['rod_constraints_test.js']
  },

  'strand_test': {
    name: 'strand_test',
    description: 'Strand Test',
    sources: ['src/strand_test.cpp'],
    useLogger: true,
    exportName: 'StrandTestModule',
    outputFiles: ['strand_test.js']
  },

  'unit_tests': {
    name: 'unit_tests',
    description: 'Gaudi Unit Tests',
    sources: ['src/unit_tests.cpp'],
    useLogger: true,
    exportName: 'UnitTests',
    outputFiles: ['unit_tests.js']
  },
  /*
    'path_test': {
    name: 'path_test',
    description: 'Path Constraint Test',
    sources: ['src/path_test.cpp'],
    useLogger: true,
    exportName: 'PathTestModule',
    outputFiles: ['path_test.js']
  },
  'examples': {
    name: 'examples',
    description: 'Hello World Example',
    sources: ['src/hello_world.cpp'],
    useLogger: false,
    exportName: 'HelloWorldModule',
    outputFiles: ['hello_world.js']
  },
  'bar_demo': {
    name: 'bar_demo',
    description: 'Bar Demo',
    sources: ['src/bar_demo.cpp'],
    useLogger: false,
    exportName: 'BarDemoModule',
    outputFiles: ['bar_demo.js']
  },
  'foo_demo': {
    name: 'foo_demo',
    description: 'Foo Demo',
    sources: ['src/foo_demo.cpp'],
    useLogger: false,
    exportName: 'FooDemoModule',
    outputFiles: ['foo_demo.js']
  },
  'test_module': {
    name: 'test_module',
    description: 'Test Module',
    sources: ['src/test_module.cpp'],
    useLogger: false,
    exportName: 'TestModule',
    outputFiles: ['test_module.js']
  }*/
};

export function getAllModules(): ModuleConfig[] {
  return Object.values(MODULE_CONFIGS);
}

export function getModuleNames(): string[] {
  return Object.keys(MODULE_CONFIGS);
} 