import React, { useState } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { Button } from './ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { Play, Pause, RotateCcw, Plus, Trash2, Activity, Clock } from 'lucide-react';
import { useWasmModule, WasmModule } from '../utils/wasmLoader';
import { GaudiLoggerRenderer, WasmLoggerAPI, WasmModule as LoggerWasmModule } from './GaudiLoggerRenderer';

interface GaudiLoggerTestModule extends WasmModule, LoggerWasmModule {
  GaudiLoggerTest: new () => GaudiLoggerTestInstance;
}

interface GaudiLoggerTestInstance extends WasmLoggerAPI {
  // Animation control
  step(): void;
  reset(): void;
  start(): void;
  stop(): void;
  get_frame_count(): number;
  get_time(): number;
  get_is_animating(): boolean;
  
  // Logger methods for testing
  add_line(x0: number, y0: number, z0: number, x1: number, y1: number, z1: number, 
           r: number, g: number, b: number, a: number): void;
  add_point(x: number, y: number, z: number, r: number, g: number, b: number, a: number): void;
  clear_logger(): void;
}

export function GaudiLoggerTest() {
  const { loading, error, module, instance } = useWasmModule<GaudiLoggerTestModule, GaudiLoggerTestInstance>(
    'gaudi_logger_test',
    (module: GaudiLoggerTestModule) => new module.GaudiLoggerTest()
  );
  
  const [isAnimating, setIsAnimating] = useState(false);
  const [frameCount, setFrameCount] = useState(0);
  const [time, setTime] = useState(0);
  const [showLines, setShowLines] = useState(true);
  const [showPoints, setShowPoints] = useState(true);

  const handleStart = () => {
    if (instance) {
      instance.start();
      setIsAnimating(true);
    }
  };

  const handleStop = () => {
    if (instance) {
      instance.stop();
      setIsAnimating(false);
    }
  };

  const handleReset = () => {
    if (instance) {
      instance.reset();
      setIsAnimating(false);
      setFrameCount(0);
      setTime(0);
    }
  };

  const handleStepOnce = () => {
    if (instance) {
      instance.step();
      setFrameCount(instance.get_frame_count());
      setTime(instance.get_time());
    }
  };

  const handleAddTestLine = () => {
    if (instance) {
      instance.add_line(0, 0, 0, 1, 1, 1, 1, 0, 0, 1);
    }
  };

  const handleAddTestPoint = () => {
    if (instance) {
      instance.add_point(0.5, 0.5, 0.5, 0, 1, 0, 1);
    }
  };

  const handleClear = () => {
    if (instance) {
      instance.clear_logger();
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading WASM module...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center">
          <p className="text-red-600 mb-4">Failed to load WASM module</p>
          <p className="text-sm text-gray-600">{error}</p>
        </div>
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <Card className="mb-4">
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <Activity className="h-5 w-5" />
            Gaudi Logger Test
          </CardTitle>
          <CardDescription>
            Test the gaudi logger system with animated visualization
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-2 mb-4">
            <Button 
              onClick={isAnimating ? handleStop : handleStart}
              variant={isAnimating ? "destructive" : "default"}
              size="sm"
            >
              {isAnimating ? <Pause className="h-4 w-4 mr-1" /> : <Play className="h-4 w-4 mr-1" />}
              {isAnimating ? 'Stop' : 'Start'}
            </Button>
            
            <Button onClick={handleReset} variant="outline" size="sm">
              <RotateCcw className="h-4 w-4 mr-1" />
              Reset
            </Button>
            
            <Button onClick={handleStepOnce} variant="outline" size="sm">
              Step
            </Button>
          </div>
          
          <div className="flex flex-wrap gap-2 mb-4">
            <Button onClick={handleAddTestLine} variant="outline" size="sm">
              <Plus className="h-4 w-4 mr-1" />
              Add Line
            </Button>
            
            <Button onClick={handleAddTestPoint} variant="outline" size="sm">
              <Plus className="h-4 w-4 mr-1" />
              Add Point
            </Button>
            
            <Button onClick={handleClear} variant="outline" size="sm">
              <Trash2 className="h-4 w-4 mr-1" />
              Clear
            </Button>
          </div>
          
          <div className="flex flex-wrap gap-2 mb-4">
            <Button 
              onClick={() => setShowLines(!showLines)}
              variant={showLines ? "default" : "outline"}
              size="sm"
            >
              {showLines ? "Hide" : "Show"} Lines
            </Button>
            
            <Button 
              onClick={() => setShowPoints(!showPoints)}
              variant={showPoints ? "default" : "outline"}
              size="sm"
            >
              {showPoints ? "Hide" : "Show"} Points
            </Button>
          </div>
          
          <div className="flex gap-4 text-sm">
            <div className="flex items-center gap-1">
              <Clock className="h-4 w-4" />
              <span>Frame: {frameCount}</span>
            </div>
            <div className="flex items-center gap-1">
              <Clock className="h-4 w-4" />
              <span>Time: {time.toFixed(2)}s</span>
            </div>
            <Badge variant={isAnimating ? "default" : "secondary"}>
              {isAnimating ? "Running" : "Stopped"}
            </Badge>
          </div>
        </CardContent>
      </Card>

      {/* 3D Scene */}
      <div className="flex-1 relative min-h-[600px]">
        <Canvas 
          camera={{ position: [3, 3, 3], fov: 60, near: 0.01, far: 1000 }}
          className="!absolute !inset-0"
          gl={{ preserveDrawingBuffer: true }}
        >
          <ambientLight intensity={0.4} />
          <pointLight position={[10, 10, 10]} />
          
          <GaudiLoggerRenderer 
            wasmInstance={instance}
            wasmModule={module}
            isPlaying={isAnimating}
            showLines={showLines}
            showPoints={showPoints}
            lineWidth={2}
            pointSize={0.1}
          />
          
          <OrbitControls />
          <gridHelper args={[10, 10]} />
          <axesHelper args={[2]} />
        </Canvas>
      </div>
    </div>
  );
}
