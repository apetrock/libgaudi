/**
 * Gaudi WebAssembly Memory Helper Utilities
 * 
 * Helper functions for working with shared memory buffers between
 * C++ and JavaScript for efficient line data transfer.
 */

export class GaudiMemoryHelper {
    constructor(wasmModule) {
        this.module = wasmModule;
    }

    /**
     * Get Float32Array view of the points buffer from C++
     * Each line consists of 2 points (start, end), each point has 3 floats (x, y, z)
     * So buffer layout is: [x0, y0, z0, x1, y1, z1, x0, y0, z0, x1, y1, z1, ...]
     */
    getPointsArray() {
        const bufferPtr = this.module.get_points_buffer();
        const bufferSize = this.module.get_points_buffer_size();
        
        if (bufferPtr === 0 || bufferSize === 0) {
            return new Float32Array(0);
        }
        
        // Create a view into the WebAssembly memory
        return new Float32Array(this.module.HEAPF32.buffer, bufferPtr, bufferSize);
    }

    /**
     * Get Float32Array view of the colors buffer from C++
     * Each line has one color with 4 floats (r, g, b, a)
     * So buffer layout is: [r0, g0, b0, a0, r1, g1, b1, a1, ...]
     */
    getColorsArray() {
        const bufferPtr = this.module.get_colors_buffer();
        const bufferSize = this.module.get_colors_buffer_size();
        
        if (bufferPtr === 0 || bufferSize === 0) {
            return new Float32Array(0);
        }
        
        return new Float32Array(this.module.HEAPF32.buffer, bufferPtr, bufferSize);
    }

    /**
     * Convert the C++ line data to Three.js BufferGeometry format
     */
    createThreeJsGeometry() {
        const points = this.getPointsArray();
        const colors = this.getColorsArray();
        const lineCount = this.module.get_line_count();

        if (lineCount === 0) {
            return null;
        }

        // Three.js expects positions and colors for each vertex
        const positions = new Float32Array(points.length);
        const vertexColors = new Float32Array(points.length / 3 * 4); // 4 components per vertex (RGBA)

        // Copy positions directly
        positions.set(points);

        // Expand colors: each line has one color, but we need color per vertex
        // Each line has 2 vertices, so duplicate each color
        for (let i = 0; i < lineCount; i++) {
            const colorIndex = i * 4;
            const vertexIndex1 = i * 8;     // First vertex of line i
            const vertexIndex2 = i * 8 + 4; // Second vertex of line i

            // Copy color to both vertices of the line
            vertexColors[vertexIndex1 + 0] = colors[colorIndex + 0]; // r
            vertexColors[vertexIndex1 + 1] = colors[colorIndex + 1]; // g
            vertexColors[vertexIndex1 + 2] = colors[colorIndex + 2]; // b
            vertexColors[vertexIndex1 + 3] = colors[colorIndex + 3]; // a

            vertexColors[vertexIndex2 + 0] = colors[colorIndex + 0]; // r
            vertexColors[vertexIndex2 + 1] = colors[colorIndex + 1]; // g
            vertexColors[vertexIndex2 + 2] = colors[colorIndex + 2]; // b
            vertexColors[vertexIndex2 + 3] = colors[colorIndex + 3]; // a
        }

        return {
            positions,
            colors: vertexColors,
            lineCount
        };
    }

    /**
     * Get current line statistics
     */
    getStats() {
        return {
            lineCount: this.module.get_line_count(),
            pointsBufferSize: this.module.get_points_buffer_size(),
            colorsBufferSize: this.module.get_colors_buffer_size(),
            totalFloats: this.module.get_points_buffer_size() + this.module.get_colors_buffer_size()
        };
    }
}
