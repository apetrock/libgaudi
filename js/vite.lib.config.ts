import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import dts from 'vite-plugin-dts';
import { resolve } from 'path';

export default defineConfig({
  plugins: [
    react(),
    dts({
      insertTypesEntry: true,
      rollupTypes: true,
    }),
  ],
  build: {
    lib: {
      entry: {
        index: resolve(__dirname, 'src/index.ts'),
        components: resolve(__dirname, 'src/components/index.ts'),
        utils: resolve(__dirname, 'src/utils/index.ts'),
      },
      formats: ['es', 'cjs'],
      fileName: (format, entryName) => `${entryName}.${format === 'es' ? 'esm' : 'js'}`,
    },
    rollupOptions: {
      external: [
        'react',
        'react-dom',
        '@react-three/fiber',
        '@react-three/drei',
        'three',
        'zustand',
        'lucide-react',
        'clsx',
        'tailwind-merge',
        'class-variance-authority',
        '@radix-ui/react-slot',
      ],
      output: {
        globals: {
          react: 'React',
          'react-dom': 'ReactDOM',
          '@react-three/fiber': 'ReactThreeFiber',
          '@react-three/drei': 'ReactThreeDrei',
          three: 'THREE',
          zustand: 'zustand',
          'lucide-react': 'lucideReact',
          clsx: 'clsx',
          'tailwind-merge': 'tailwindMerge',
          'class-variance-authority': 'classVarianceAuthority',
          '@radix-ui/react-slot': 'radixSlot',
        },
      },
    },
    outDir: 'dist',
    sourcemap: true,
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, 'src'),
    },
  },
}); 