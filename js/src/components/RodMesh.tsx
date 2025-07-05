import React, { useMemo, useRef, useEffect } from 'react';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';
import type { RodSimulationData, RodMaterialProps } from '../types/simulation';

interface RodMeshProps extends RodMaterialProps {
  data: RodSimulationData | null;
  tubularSegments?: number;
  radialSegments?: number;
  radius?: number;
}

/**
 * Component that renders the dynamic rod as a tube geometry
 * Updates vertex positions and colors based on simulation data
 */
export function RodMesh({
  data,
  tubularSegments = 256,
  radialSegments = 8,
  radius = 0.02,
  baseColor = '#4a90e2',
  growthColorLow = '#2ecc71',
  growthColorHigh = '#e74c3c',
  metalness = 0.1,
  roughness = 0.4,
  wireframe = false
}: RodMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null);
  const geometryRef = useRef<THREE.TubeGeometry | null>(null);
  const materialRef = useRef<THREE.MeshStandardMaterial>(null);

  // Create curve from rod vertices
  const curve = useMemo(() => {
    if (!data || data.vertices.length === 0) {
      // Default spiral curve for initial state
      return new THREE.CatmullRomCurve3([
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(1, 0, 0),
        new THREE.Vector3(1, 1, 0),
        new THREE.Vector3(0, 1, 0)
      ], true);
    }

    const points: THREE.Vector3[] = [];
    for (let i = 0; i < data.vertices.length; i += 3) {
      points.push(new THREE.Vector3(
        data.vertices[i],
        data.vertices[i + 1], 
        data.vertices[i + 2]
      ));
    }

    return new THREE.CatmullRomCurve3(points, true);
  }, [data]);

  // Create geometry
  const geometry = useMemo(() => {
    const geom = new THREE.TubeGeometry(curve, tubularSegments, radius, radialSegments, true);
    geometryRef.current = geom;
    return geom;
  }, [curve, tubularSegments, radius, radialSegments]);

  // Create material with vertex colors for growth visualization
  const material = useMemo(() => {
    const mat = new THREE.MeshStandardMaterial({
      color: baseColor,
      metalness,
      roughness,
      wireframe,
      vertexColors: true
    });
    materialRef.current = mat;
    return mat;
  }, [baseColor, metalness, roughness, wireframe]);

  // Update vertex colors based on growth weights
  useEffect(() => {
    if (!data || !geometryRef.current) return;

    const geometry = geometryRef.current;
    const colorAttribute = geometry.getAttribute('color') as THREE.BufferAttribute;
    
    if (!colorAttribute) {
      // Create color attribute if it doesn't exist
      const colors = new Float32Array(geometry.attributes.position.count * 3);
      geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    }

    const colors = geometry.getAttribute('color') as THREE.BufferAttribute;
    const lowColor = new THREE.Color(growthColorLow);
    const highColor = new THREE.Color(growthColorHigh);
    const tempColor = new THREE.Color();

    // Map growth weights to colors
    for (let i = 0; i < colors.count; i++) {
      // Get the corresponding growth weight (accounting for tube geometry structure)
      const rodVertexIndex = Math.floor(i / radialSegments) % data.growthWeights.length;
      const weight = data.growthWeights[rodVertexIndex] || 0;
      
      // Interpolate between low and high growth colors
      const normalizedWeight = THREE.MathUtils.clamp(weight / 4.0, 0, 1);
      tempColor.lerpColors(lowColor, highColor, normalizedWeight);
      
      colors.setXYZ(i, tempColor.r, tempColor.g, tempColor.b);
    }

    colors.needsUpdate = true;
  }, [data, growthColorLow, growthColorHigh, radialSegments]);

  // Update geometry when curve changes
  useFrame(() => {
    if (!data || !geometryRef.current) return;

    // Update the curve and regenerate geometry if needed
    // This is a performance consideration - we might want to update positions directly
    // instead of recreating the entire geometry each frame
  });

  return (
    <mesh ref={meshRef} geometry={geometry} material={material}>
      {/* Optional: Add custom shader material for more advanced effects */}
    </mesh>
  );
}

// Helper component for SDF visualization
export function SdfVisualization({ 
  frame, 
  opacity = 0.3, 
  visible = true 
}: {
  frame: number;
  opacity?: number;
  visible?: boolean;
}) {
  const sphere1Ref = useRef<THREE.Mesh>(null);
  const sphere2Ref = useRef<THREE.Mesh>(null);

  // Determine which SDF to show based on frame
  const showMultiSphere = Math.floor(frame / 400) % 2 === 1;

  // Fibonacci sphere positions (matching C++ implementation)
  const fibPositions = useMemo(() => {
    const positions: THREE.Vector3[] = [];
    const golden = 0.5 * (1.0 + Math.sqrt(5));
    const r = 1.5;
    const N = 13;

    for (let i = 0; i < N; i++) {
      const theta = 2.0 * Math.PI * i / golden;
      const phi = Math.acos(1.0 - 2.0 * (i + 0.5) / N);
      const x = r * Math.cos(theta) * Math.sin(phi);
      const y = r * Math.sin(theta) * Math.sin(phi);
      const z = r * Math.cos(phi);
      positions.push(new THREE.Vector3(x, y, z));
    }

    return positions;
  }, []);

  return (
    <group visible={visible}>
      {/* Single sphere SDF */}
      <mesh ref={sphere1Ref} visible={!showMultiSphere}>
        <sphereGeometry args={[1.26, 32, 32]} />
        <meshStandardMaterial 
          color="#3498db" 
          transparent 
          opacity={opacity}
          wireframe
        />
      </mesh>

      {/* Multi-sphere SDF */}
      <group visible={showMultiSphere}>
        {fibPositions.map((pos, index) => (
          <mesh key={index} position={[pos.x, pos.y, pos.z]}>
            <sphereGeometry args={[0.5, 16, 16]} />
            <meshStandardMaterial 
              color="#e67e22" 
              transparent 
              opacity={opacity}
              wireframe
            />
          </mesh>
        ))}
      </group>
    </group>
  );
}
