import React, { useState, useRef, useCallback } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { Button } from './ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { Play, Pause, RotateCcw, Activity, Clock, Eye, EyeOff, Video, Download } from 'lucide-react';
import { useWasmModule, WasmModule } from '../utils/wasmLoader';
import { GaudiLoggerRenderer, WasmLoggerAPI, WasmModule as LoggerWasmModule } from './GaudiLoggerRenderer';
import { VideoRecorder } from './VideoRecorder';

interface StrandTestModule extends WasmModule, LoggerWasmModule {
  StrandTest: new () => StrandTestInstance;
}

interface StrandTestInstance extends WasmLoggerAPI {
  // Animation control
  step(): void;
  reset(): void;
  get_frame_count(): number;
  get_time(): number;
}

// Component to handle simulation stepping
function SimulationController({ 
  wasmInstance, 
  isPlaying, 
  onFrameUpdate 
}: { 
  wasmInstance: StrandTestInstance | null;
  isPlaying: boolean;
  onFrameUpdate: (frame: number) => void;
}) {
  useFrame(() => {
    if (!wasmInstance || !isPlaying) return;
    
    // Step the simulation
    wasmInstance.step();
    
    // Update frame count
    onFrameUpdate(wasmInstance.get_frame_count());
  });

  return null; // This component doesn't render anything
}

export function StrandTest() {
  const instanceFactory = useCallback((module: StrandTestModule) => {
    return new module.StrandTest();
  }, []);

  const { loading, error, module, instance } = useWasmModule<StrandTestModule, StrandTestInstance>(
    'strand_test',
    instanceFactory
  );
  
  const [isPlaying, setIsPlaying] = useState(false);
  const [frame, setFrame] = useState(0);
  const [showDebugLines, setShowDebugLines] = useState(true);
  const [showDebugPoints, setShowDebugPoints] = useState(true);
  const [recordingBlob, setRecordingBlob] = useState<Blob | null>(null);
  
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const handlePlayPause = () => {
    setIsPlaying(!isPlaying);
  };

  const handleReset = () => {
    if (instance) {
      instance.reset();
      setFrame(0);
    }
  };

  const handleStepOnce = () => {
    if (instance) {
      instance.step();
      setFrame(instance.get_frame_count());
    }
  };

  const handleRecordingComplete = (blob: Blob) => {
    setRecordingBlob(blob);
  };

  const handleDownloadRecording = () => {
    if (recordingBlob) {
      const url = URL.createObjectURL(recordingBlob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `rod-strand-test-${Date.now()}.webm`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
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
          <p className="text-xs text-gray-500 mt-2">
            Build the module with: cd wasm/rod_strand && ./build.bat
          </p>
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
            Rod Strand Test
          </CardTitle>
          <CardDescription>
            Test rod strand dynamics with logger visualization and video recording
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-2 mb-4">
            <Button 
              onClick={handlePlayPause}
              variant={isPlaying ? "destructive" : "default"}
              size="sm"
            >
              {isPlaying ? <Pause className="h-4 w-4 mr-1" /> : <Play className="h-4 w-4 mr-1" />}
              {isPlaying ? 'Pause' : 'Play'}
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
            <Button 
              onClick={() => setShowDebugLines(!showDebugLines)}
              variant={showDebugLines ? "default" : "outline"}
              size="sm"
            >
              {showDebugLines ? <Eye className="h-4 w-4 mr-1" /> : <EyeOff className="h-4 w-4 mr-1" />}
              Lines
            </Button>
            
            <Button 
              onClick={() => setShowDebugPoints(!showDebugPoints)}
              variant={showDebugPoints ? "default" : "outline"}
              size="sm"
            >
              {showDebugPoints ? <Eye className="h-4 w-4 mr-1" /> : <EyeOff className="h-4 w-4 mr-1" />}
              Points
            </Button>
          </div>

          {/* Video Recording Controls */}
          <div className="flex flex-wrap gap-2 mb-4">
            <VideoRecorder 
              targetRef={canvasRef}
              quality="high"
              framerate={60}
              onRecordingComplete={handleRecordingComplete}
            >
              {({ startRecording, stopRecording, isRecording }) => (
                <>
                  <Button 
                    onClick={isRecording ? stopRecording : startRecording}
                    variant={isRecording ? "destructive" : "outline"}
                    size="sm"
                  >
                    <Video className="h-4 w-4 mr-1" />
                    {isRecording ? 'Stop Recording' : 'Start Recording'}
                  </Button>
                </>
              )}
            </VideoRecorder>
            
            {recordingBlob && (
              <Button onClick={handleDownloadRecording} variant="outline" size="sm">
                <Download className="h-4 w-4 mr-1" />
                Download Recording
              </Button>
            )}
          </div>
          
          <div className="flex gap-4 text-sm">
            <div className="flex items-center gap-1">
              <Clock className="h-4 w-4" />
              <span>Frame: {frame}</span>
            </div>
            <div className="flex items-center gap-1">
              <Clock className="h-4 w-4" />
              <span>Time: {instance?.get_time().toFixed(2)}s</span>
            </div>
            <Badge variant={isPlaying ? "default" : "secondary"}>
              {isPlaying ? "Running" : "Stopped"}
            </Badge>
          </div>
        </CardContent>
      </Card>

      {/* 3D Scene */}
      <div className="flex-1 relative min-h-[600px]">
        <Canvas 
          ref={canvasRef}
          camera={{ position: [3, 3, 3], fov: 60, near: 0.01, far: 1000 }}
          className="!absolute !inset-0"
          gl={{ preserveDrawingBuffer: true }}
        >
          <ambientLight intensity={0.4} />
          <pointLight position={[10, 10, 10]} />
          
          <SimulationController 
            wasmInstance={instance}
            isPlaying={isPlaying}
            onFrameUpdate={setFrame}
          />
          
          <GaudiLoggerRenderer 
            wasmInstance={instance}
            wasmModule={module}
            isPlaying={isPlaying}
            showLines={showDebugLines}
            showPoints={showDebugPoints}
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