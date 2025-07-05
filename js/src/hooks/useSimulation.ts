import { useState, useCallback, useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import type { 
  SimulationParams, 
  RodSimulationData, 
  WasmSimulation 
} from '../types/simulation';

const DEFAULT_PARAMS: SimulationParams = {
  growthRate: 1.0,
  constraintStrength: 0.1,
  collisionStrength: 1.0,
  timeStep: 0.05,
  sdfTransitionFrame: 400,
  bendingStiffness: 0.05,
  stretchStiffness: 0.06,
};

interface UseSimulationReturn {
  data: RodSimulationData | null;
  params: SimulationParams;
  isPlaying: boolean;
  frame: number;
  play: () => void;
  pause: () => void;
  reset: () => void;
  step: () => void;
  updateParams: (newParams: Partial<SimulationParams>) => void;
}

/**
 * Hook to manage simulation state and animation loop
 */
export function useSimulation(
  wasmSimulation: WasmSimulation | null,
  initialParams?: Partial<SimulationParams>,
  onFrameUpdate?: (data: RodSimulationData) => void
): UseSimulationReturn {
  const [params, setParams] = useState<SimulationParams>(() => ({
    ...DEFAULT_PARAMS,
    ...initialParams
  }));
  
  const [isPlaying, setIsPlaying] = useState(false);
  const [frame, setFrame] = useState(0);
  const [data, setData] = useState<RodSimulationData | null>(null);
  
  const frameRef = useRef(0);
  const paramsRef = useRef(params);
  
  // Keep params ref in sync
  paramsRef.current = params;
  // Extract data from WASM simulation
  const extractSimulationData = useCallback((simulation: WasmSimulation, currentFrame: number): RodSimulationData => {
    const verticesArray = simulation.getVertices();
    const normalsArray = simulation.getNormals();
    const tangentsArray = simulation.getTangets(); // Note: matches our C++ typo
    const growthWeightsArray = simulation.getGrowthWeights();
    
    // Convert JS arrays to Float32Arrays for consistent interface
    const vertices = new Float32Array(verticesArray);
    const normals = new Float32Array(normalsArray);
    const tangents = new Float32Array(tangentsArray);
    const growthWeights = new Float32Array(growthWeightsArray);
    
    const vertexCount = simulation.getVertexCount();
    const totalLength = simulation.getTotalLength();

    return {
      vertices,
      normals,
      tangents,
      growthWeights,
      frame: currentFrame,
      totalLength,
      vertexCount
    };
  }, []);

  // Animation loop using R3F's useFrame
  useFrame(() => {
    if (!wasmSimulation || !isPlaying) return;

    try {
      // Step the simulation
      wasmSimulation.step(frameRef.current);
      
      // Extract current data
      const currentData = extractSimulationData(wasmSimulation, frameRef.current);
      setData(currentData);
      setFrame(frameRef.current);
      
      // Notify parent component
      if (onFrameUpdate) {
        onFrameUpdate(currentData);
      }
      
      frameRef.current++;
    } catch (error) {
      console.error('Simulation step failed:', error);
      setIsPlaying(false);
    }
  });

  // Initialize simulation data when WASM simulation becomes available
  useState(() => {
    if (wasmSimulation && !data) {
      try {
        const initialData = extractSimulationData(wasmSimulation, 0);
        setData(initialData);
      } catch (error) {
        console.error('Failed to extract initial simulation data:', error);
      }
    }
  });

  const play = useCallback(() => {
    setIsPlaying(true);
  }, []);

  const pause = useCallback(() => {
    setIsPlaying(false);
  }, []);

  const reset = useCallback(() => {
    if (!wasmSimulation) return;
    
    setIsPlaying(false);
    frameRef.current = 0;
    setFrame(0);
    
    try {
      wasmSimulation.reset();
      const resetData = extractSimulationData(wasmSimulation, 0);
      setData(resetData);
      
      if (onFrameUpdate) {
        onFrameUpdate(resetData);
      }
    } catch (error) {
      console.error('Failed to reset simulation:', error);
    }
  }, [wasmSimulation, extractSimulationData, onFrameUpdate]);

  const step = useCallback(() => {
    if (!wasmSimulation) return;
    
    try {
      wasmSimulation.step(frameRef.current);
      const stepData = extractSimulationData(wasmSimulation, frameRef.current);
      setData(stepData);
      setFrame(frameRef.current);
      
      if (onFrameUpdate) {
        onFrameUpdate(stepData);
      }
      
      frameRef.current++;
    } catch (error) {
      console.error('Manual simulation step failed:', error);
    }
  }, [wasmSimulation, extractSimulationData, onFrameUpdate]);
  const updateParams = useCallback((newParams: Partial<SimulationParams>) => {
    const updatedParams = { ...params, ...newParams };
    setParams(updatedParams);
    
    // Update WASM simulation parameters - map to our specific setters
    if (wasmSimulation) {
      try {
        if (newParams.growthRate !== undefined) {
          wasmSimulation.setGrowthRate(newParams.growthRate);
        }
        if (newParams.constraintStrength !== undefined) {
          wasmSimulation.setConstraintStrength(newParams.constraintStrength);
        }
        // TODO: Add other parameter setters as we expand the C++ interface
      } catch (error) {
        console.error('Failed to update simulation parameters:', error);
      }
    }
  }, [params, wasmSimulation]);

  return {
    data,
    params,
    isPlaying,
    frame,
    play,
    pause,
    reset,
    step,
    updateParams
  };
}
