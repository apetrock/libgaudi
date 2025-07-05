import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import dts from 'vite-plugin-dts'
import path from 'path'
import { fileURLToPath } from 'url'
import { dirname } from 'path'

const __filename = fileURLToPath(import.meta.url)
const __dirname = dirname(__filename)

// Custom plugin to handle WASM files properly
const wasmPlugin = () => {
  return {
    name: 'wasm-mime-type',
    configureServer(server) {
      server.middlewares.use((req, res, next) => {
        if (req.url?.endsWith('.wasm')) {
          res.setHeader('Content-Type', 'application/wasm')
        }
        next()
      })
    }
  }
}

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  if (mode === 'lib') {
    // Library build configuration
    return {
      plugins: [
        react(),
        dts({
          insertTypesEntry: true,
        })
      ],
      build: {
        lib: {
          entry: path.resolve(__dirname, 'src/index.ts'),
          name: 'RodSimulationComponent',
          formats: ['es', 'umd'],
          fileName: (format) => `index.${format === 'es' ? 'esm' : format}.js`
        },
        rollupOptions: {
          external: [
            'react',
            'react-dom',
            '@react-three/fiber',
            '@react-three/drei',
            'three'
          ],
          output: {
            globals: {
              react: 'React',
              'react-dom': 'ReactDOM',
              '@react-three/fiber': 'ReactThreeFiber',
              '@react-three/drei': 'ReactThreeDrei',
              three: 'THREE'
            }
          }
        }
      }
    }
  }
  
  // Development configuration - simplified like Henson
  return {
    plugins: [
      react(),
      wasmPlugin()
    ],
    resolve: {
      alias: {
        '@': path.resolve(__dirname, './src'),
      },
    },
    server: {
      port: 5173,
      headers: {
        'Cross-Origin-Embedder-Policy': 'require-corp',
        'Cross-Origin-Opener-Policy': 'same-origin',
      },
      fs: {
        allow: ['..']
      }
    },
    appType: 'spa',
    assetsInclude: ['**/*.wasm'],
    optimizeDeps: {
      exclude: ['*.wasm']
    },
    // Ensure public directory is served correctly
    publicDir: 'public'
  }
})
