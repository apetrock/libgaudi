import React, { useState, useEffect, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { Button } from './ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { Play, Pause, RotateCcw, Plus, Trash2, Activity, Clock } from 'lucide-react';
import { loadWasmModule } from '../utils/wasmLoader';

interface GaudiLoggerTestModule {
  GaudiLoggerTest: new () => GaudiLoggerTestInstance;
  HEAPF64: Float64Array;
}

interface GaudiLoggerTestInstance {
  // Animation control
  step(): void;
  reset(): void;
  get_frame_count(): number;
  get_time(): number;
  
  // Line data access
  get_line_count(): number;
  get_line_positions_ptr(): number;
  get_line_positions_size(): number;
  get_line_colors_ptr(): number;
  get_line_colors_size(): number;
  
  // Point data access
  get_point_count(): number;
  get_point_positions_ptr(): number;
  get_point_positions_size(): number;
  get_point_colors_ptr(): number;
  get_point_colors_size(): number;
  
  // Manual testing
  add_test_line(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number,
                r: number, g: number, b: number, a: number): void;
  add_test_point(x: number, y: number, z: number, r: number, g: number, b: number, a: number): void;
  clear_logger(): void;
}

interface WasmState {
  loading: boolean;
  error: string | null;
  module: GaudiLoggerTestModule | null;
  instance: GaudiLoggerTestInstance | null;
}

function GaudiLoggerGeometry({ 
  wasmModule, 
  wasmInstance, 
  isAnimating 
}: { 
  wasmModule: GaudiLoggerTestModule | null;
  wasmInstance: GaudiLoggerTestInstance | null;
  isAnimating: boolean;
}) {
  const linesRef = useRef<THREE.LineSegments>(null);
  const pointsRef = useRef<THREE.Points>(null);

  useFrame(() => {
    if (!wasmModule || !wasmInstance) return;

    // Step the animation
    if (isAnimating) {
      wasmInstance.step();
    }

    // Update Lines
    const lineCount = wasmInstance.get_line_count();
    if (linesRef.current) {
      if (lineCount > 0) {
        const positionsPtr = wasmInstance.get_line_positions_ptr();
        const colorsPtr = wasmInstance.get_line_colors_ptr();
        
        // Get line data from WASM heap
        const positions = new Float64Array(wasmModule.HEAPF64.buffer, positionsPtr, lineCount * 2 * 3);
        const colors = new Float64Array(wasmModule.HEAPF64.buffer, colorsPtr, lineCount * 2 * 4);
        
        // Convert to Three.js format
        const threePositions = new Float32Array(positions);
        const threeColors = new Float32Array(lineCount * 2 * 3);

        for (let i = 0; i < lineCount * 2; i++) {
          threeColors[i * 3 + 0] = colors[i * 4 + 0];
          threeColors[i * 3 + 1] = colors[i * 4 + 1];
          threeColors[i * 3 + 2] = colors[i * 4 + 2];
        }

        linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(threePositions, 3));
        linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(threeColors, 3));
        linesRef.current.geometry.attributes.position.needsUpdate = true;
        linesRef.current.geometry.attributes.color.needsUpdate = true;
      } else {
        linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 3));
      }
    }

    // Update Points
    const pointCount = wasmInstance.get_point_count();
    if (pointsRef.current) {
      if (pointCount > 0) {
        const positionsPtr = wasmInstance.get_point_positions_ptr();
        const colorsPtr = wasmInstance.get_point_colors_ptr();
        
        const positions = new Float64Array(wasmModule.HEAPF64.buffer, positionsPtr, pointCount * 3);
        const colors = new Float64Array(wasmModule.HEAPF64.buffer, colorsPtr, pointCount * 4);

        const threePositions = new Float32Array(positions);
        const threeColors = new Float32Array(pointCount * 3);

        for (let i = 0; i < pointCount; i++) {
          threeColors[i * 3 + 0] = colors[i * 4 + 0];
          threeColors[i * 3 + 1] = colors[i * 4 + 1];
          threeColors[i * 3 + 2] = colors[i * 4 + 2];
        }

        pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(threePositions, 3));
        pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(threeColors, 3));
        pointsRef.current.geometry.attributes.position.needsUpdate = true;
        pointsRef.current.geometry.attributes.color.needsUpdate = true;
      } else {
        pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 3));
      }
    }
  });

  return (
    <>
      <lineSegments ref={linesRef}>
        <bufferGeometry />
        <lineBasicMaterial vertexColors={true} />
      </lineSegments>
      <points ref={pointsRef}>
        <bufferGeometry />
        <pointsMaterial size={5} vertexColors={true} />
      </points>
    </>
  );
}

function ControlPanel({ 
  wasmInstance, 
  isAnimating, 
  setIsAnimating 
}: { 
  wasmInstance: GaudiLoggerTestInstance | null;
  isAnimating: boolean;
  setIsAnimating: (value: boolean) => void;
}) {
  const [frameCount, setFrameCount] = useState(0);
  const [time, setTime] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => {
      if (wasmInstance) {
        setFrameCount(wasmInstance.get_frame_count());
        setTime(wasmInstance.get_time());
      }
    }, 100);

    return () => clearInterval(interval);
  }, [wasmInstance]);

  const handleReset = () => {
    if (wasmInstance) {
      wasmInstance.reset();
      setFrameCount(0);
      setTime(0);
    }
  };

  const handleAddTestLine = () => {
    if (wasmInstance) {
      const x1 = (Math.random() - 0.5) * 4;
      const y1 = (Math.random() - 0.5) * 4;
      const z1 = (Math.random() - 0.5) * 4;
      const x2 = x1 + (Math.random() - 0.5) * 2;
      const y2 = y1 + (Math.random() - 0.5) * 2;
      const z2 = z1 + (Math.random() - 0.5) * 2;
      wasmInstance.add_test_line(x1, y1, z1, x2, y2, z2, Math.random(), Math.random(), Math.random(), 1.0);
    }
  };

  const handleClearLogger = () => {
    if (wasmInstance) {
      wasmInstance.clear_logger();
    }
  };

  return (
    <Card className="absolute top-4 left-4 w-80 bg-background/80 backdrop-blur-sm border-border/50">
      <CardHeader className="pb-3">
        <CardTitle className="flex items-center space-x-2 text-lg">
          <Activity className="h-5 w-5 text-primary" />
          <span>Gaudi Logger Test</span>
          <Badge variant="secondary" className="ml-auto">
            Phase 1.5
          </Badge>
        </CardTitle>
        <CardDescription>
          Real-time C++ logger interface testing
        </CardDescription>
      </CardHeader>
      
      <CardContent className="space-y-4">
        {/* Stats */}
        <div className="grid grid-cols-2 gap-3">
          <div className="flex items-center space-x-2 p-2 bg-muted/50 rounded-md">
            <Activity className="h-4 w-4 text-primary" />
            <div className="text-sm">
              <div className="font-medium">{frameCount}</div>
              <div className="text-muted-foreground text-xs">Frames</div>
            </div>
          </div>
          <div className="flex items-center space-x-2 p-2 bg-muted/50 rounded-md">
            <Clock className="h-4 w-4 text-primary" />
            <div className="text-sm">
              <div className="font-medium">{time.toFixed(2)}s</div>
              <div className="text-muted-foreground text-xs">Time</div>
            </div>
          </div>
        </div>
        
        {/* Animation Controls */}
        <div className="flex space-x-2">
          <Button
            onClick={() => setIsAnimating(!isAnimating)}
            variant={isAnimating ? "destructive" : "default"}
            size="sm"
            className="flex-1"
          >
            {isAnimating ? (
              <>
                <Pause className="h-4 w-4 mr-2" />
                Pause
              </>
            ) : (
              <>
                <Play className="h-4 w-4 mr-2" />
                Play
              </>
            )}
          </Button>
          
          <Button
            onClick={handleReset}
            variant="outline"
            size="sm"
          >
            <RotateCcw className="h-4 w-4" />
          </Button>
        </div>
        
        {/* Manual Controls */}
        <div className="space-y-2">
          <div className="text-sm font-medium text-muted-foreground">Manual Testing</div>
          <div className="flex space-x-2">
            <Button
              onClick={handleAddTestLine}
              variant="secondary"
              size="sm"
              className="flex-1"
            >
              <Plus className="h-4 w-4 mr-2" />
              Add Line
            </Button>
            
            <Button
              onClick={handleClearLogger}
              variant="outline"
              size="sm"
            >
              <Trash2 className="h-4 w-4" />
            </Button>
          </div>
        </div>
        
        {/* Features List */}
        <div className="pt-2 border-t border-border/50">
          <div className="text-xs text-muted-foreground space-y-1">
            <div>✨ Rotating coordinate frame</div>
            <div>🌊 Animated sine wave</div>
            <div>📦 Orbiting box animation</div>
            <div>🎛️ Manual line addition</div>
            <div>⚡ Real-time gaudi::logger interface</div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

export function GaudiLoggerTest() {
  const [wasmState, setWasmState] = useState<WasmState>({
    loading: true,
    error: null,
    module: null,
    instance: null
  });
  const [isAnimating, setIsAnimating] = useState(true);

  useEffect(() => {
    const loadWasm = async () => {
      try {
        const module = await loadWasmModule('gaudi_logger_test');
        const instance = new module.GaudiLoggerTest();
        
        setWasmState({
          loading: false,
          error: null,
          module,
          instance
        });
      } catch (error) {
        setWasmState({
          loading: false,
          error: error instanceof Error ? error.message : 'Unknown error',
          module: null,
          instance: null
        });
      }
    };

    loadWasm();
  }, []);

  if (wasmState.loading) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-background">
        <Card className="w-96">
          <CardHeader>
            <CardTitle className="flex items-center space-x-2">
              <Activity className="h-5 w-5 text-primary animate-pulse" />
              <span>Loading Gaudi Logger Test</span>
            </CardTitle>
            <CardDescription>
              Initializing WASM module...
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="w-full bg-muted rounded-full h-2">
              <div className="bg-primary h-2 rounded-full animate-pulse w-3/4"></div>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  if (wasmState.error) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-background">
        <Card className="w-96 border-destructive/50">
          <CardHeader>
            <CardTitle className="text-destructive">
              Failed to load WASM module
            </CardTitle>
            <CardDescription>
              {wasmState.error}
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="text-xs text-muted-foreground bg-muted p-3 rounded-md font-mono">
              Make sure the WASM module is built:<br />
              <code>cd wasm && build_gaudi_logger_test.bat</code>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  return (
    <div className="w-full h-full relative">
      <Canvas
        camera={{ position: [8, 6, 8], fov: 60 }}
        className="bg-background"
      >
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} />
        
        <GaudiLoggerGeometry 
          wasmModule={wasmState.module}
          wasmInstance={wasmState.instance}
          isAnimating={isAnimating}
        />
        
        <OrbitControls 
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
        />
        
        {/* Grid helper */}
        <gridHelper args={[20, 20, '#444', '#222']} />
        <axesHelper args={[3]} />
      </Canvas>
      
      <ControlPanel 
        wasmInstance={wasmState.instance}
        isAnimating={isAnimating}
        setIsAnimating={setIsAnimating}
      />
    </div>
  );
}
