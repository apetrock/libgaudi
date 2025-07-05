import { useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';

interface BoxMeshProps {
  color?: string;
  rotationSpeed?: number;
  position?: [number, number, number];
}

/**
 * Simple animated box component for validation
 */
export function BoxMesh({ 
  color = '#4a90e2', 
  rotationSpeed = 1.0,
  position = [0, 0, 0] 
}: BoxMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null);
  
  // Animation loop
  useFrame((state, delta) => {
    if (meshRef.current) {
      meshRef.current.rotation.x += delta * rotationSpeed;
      meshRef.current.rotation.y += delta * rotationSpeed * 0.7;
      
      // Add a subtle floating animation
      meshRef.current.position.y = position[1] + Math.sin(state.clock.elapsedTime) * 0.1;
    }
  });

  return (
    <mesh ref={meshRef} position={position}>
      <boxGeometry args={[1, 1, 1]} />
      <meshStandardMaterial color={color} />
    </mesh>
  );
}
