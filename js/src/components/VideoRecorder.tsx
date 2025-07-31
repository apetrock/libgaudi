import React, { useRef, useCallback, useState, useEffect } from 'react';

interface VideoRecorderProps {
  targetRef: React.RefObject<HTMLElement>;
  onRecordingChange?: (isRecording: boolean) => void;
  onRecordingComplete?: (blob: Blob) => void;
  quality?: 'low' | 'medium' | 'high';
  framerate?: number;
  codec?: 'vp8' | 'vp9' | 'h264';
  children?: React.ReactNode | ((props: { startRecording: () => void; stopRecording: () => void; isRecording: boolean }) => React.ReactNode);
}

interface RecordingConfig {
  width: number;
  height: number;
  bitrate: number;
  framerate: number;
  codec: string;
}

/**
 * Video Recorder Component using WebCodecs API
 * 
 * Usage:
 * <VideoRecorder targetRef={canvasRef}>
 *   <button onClick={() => startRecording()}>Start</button>
 *   <button onClick={() => stopRecording()}>Stop</button>
 * </VideoRecorder>
 */
export function VideoRecorder({
  targetRef,
  onRecordingChange,
  onRecordingComplete,
  quality = 'medium',
  framerate = 60,
  codec = 'vp8',
  children
}: VideoRecorderProps) {
  const [isRecording, setIsRecording] = useState(false);
  const [isSupported, setIsSupported] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const encoderRef = useRef<VideoEncoder | null>(null);
  const chunksRef = useRef<EncodedVideoChunk[]>([]);
  const animationFrameRef = useRef<number | null>(null);
  const startTimeRef = useRef<number>(0);
  const [config, setConfig] = useState<RecordingConfig | null>(null);
  // Check WebCodecs support
  useEffect(() => {
    const supported = typeof VideoEncoder !== 'undefined' && typeof VideoFrame !== 'undefined';
    setIsSupported(supported);
    if (!supported) {
      setError('WebCodecs API not supported in this browser');
    }
  }, []);

  // Quality presets
  const getQualityConfig = (quality: string): RecordingConfig => {
    const configs = {
      low: { bitrate: 500000, framerate: 30 },
      medium: { bitrate: 1000000, framerate: 60 },
      high: { bitrate: 2000000, framerate: 60 }
    };
    
    const config = configs[quality as keyof typeof configs] || configs.medium;
    return {
      width: targetRef.current?.offsetWidth || 1920,
      height: targetRef.current?.offsetHeight || 1080,
      bitrate: config.bitrate,
      framerate: config.framerate,
      codec
    };
  };


  useEffect(() => {
    
      // Start frame capture loop
      const captureFrame = () => {
        if (!isRecording || !encoderRef.current || !targetRef.current) return;

        try {
          // Create a canvas from the target element
          const canvas = document.createElement('canvas');
          const ctx = canvas.getContext('2d');
          if (!ctx) return;
          if(!config) return;
          
          canvas.width = config.width;
          canvas.height = config.height;

          // Use html2canvas or similar to capture the div
          // For now, we'll assume the target is a canvas
          if (targetRef.current instanceof HTMLCanvasElement) {
            ctx.drawImage(targetRef.current, 0, 0, config.width, config.height);
          } else {
            // For non-canvas elements, you'd need html2canvas
            // This is a placeholder - you'd need to implement div capture
            ctx.fillStyle = '#000';
            ctx.fillRect(0, 0, config.width, config.height);
          }

          const frame = new VideoFrame(canvas, {
            timestamp: performance.now() - startTimeRef.current
          });
          encoderRef.current.encode(frame);
          frame.close();

          animationFrameRef.current = requestAnimationFrame(captureFrame);
        } catch (err) {
          setError(`Frame capture error: ${err}`);
        }
      };

      captureFrame();

  }, [isRecording, config]);

  // Start recording
  const startRecording = useCallback(async () => {
    if (!isSupported || !targetRef.current) {
      setError('WebCodecs not supported or target element not found');
      return;
    }

    try {
      setError(null);
      chunksRef.current = [];
      
      const config = getQualityConfig(quality);
      
      encoderRef.current = new VideoEncoder({
        output: (chunk: EncodedVideoChunk) => {
          chunksRef.current.push(chunk);
        },
        error: (error: Error) => {
          setError(`Encoder error: ${error.message}`);
          setIsRecording(false);
        }
      });

      await encoderRef.current.configure({
        codec: config.codec as any,
        width: config.width,
        height: config.height,
        bitrate: config.bitrate,
        framerate: config.framerate
      });

      setIsRecording(true);
      setConfig(config);
      startTimeRef.current = performance.now();
      onRecordingChange?.(true);

    } catch (err) {
      setError(`Failed to start recording: ${err}`);
    }
  }, [isSupported, targetRef, quality, codec, framerate, isRecording, onRecordingChange]);

  // Stop recording
  const stopRecording = useCallback(async () => {
    if (!encoderRef.current) return;

    try {
      setIsRecording(false);
      onRecordingChange?.(false);

      if (animationFrameRef.current) {
        cancelAnimationFrame(animationFrameRef.current);
        animationFrameRef.current = null;
      }

      await encoderRef.current.flush();
      encoderRef.current.close();
      encoderRef.current = null;

      // Create blob from chunks
      const blob = new Blob(chunksRef.current as unknown as BlobPart[], { 
        type: `video/${codec === 'h264' ? 'mp4' : 'webm'}` 
      });
      
      chunksRef.current = [];
      onRecordingComplete?.(blob);

    } catch (err) {
      setError(`Failed to stop recording: ${err}`);
    }
  }, [codec, onRecordingChange, onRecordingComplete]);

  // Download helper
  const downloadRecording = useCallback((blob: Blob, filename = 'recording.webm') => {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      if (animationFrameRef.current) {
        cancelAnimationFrame(animationFrameRef.current);
      }
      if (encoderRef.current) {
        encoderRef.current.close();
      }
    };
  }, []);

  return (
    <div className="video-recorder">
      {error && (
        <div style={{
          background: '#fee',
          color: '#c33',
          padding: '8px',
          borderRadius: '4px',
          marginBottom: '8px',
          fontSize: '12px'
        }}>
          {error}
        </div>
      )}
      
      {!isSupported && (
        <div style={{
          background: '#fef',
          color: '#33c',
          padding: '8px',
          borderRadius: '4px',
          marginBottom: '8px',
          fontSize: '12px'
        }}>
          WebCodecs not supported. Try Chrome 94+ or Firefox 108+.
        </div>
      )}

      {children && typeof children === 'function' ? children({ startRecording, stopRecording, isRecording }) : children}

      {/* Recording status indicator */}
      {isRecording && (
        <div style={{
          position: 'fixed',
          top: '20px',
          right: '20px',
          background: '#f44336',
          color: 'white',
          padding: '8px 12px',
          borderRadius: '4px',
          fontSize: '12px',
          zIndex: 1000,
          display: 'flex',
          alignItems: 'center',
          gap: '8px'
        }}>
          <div style={{
            width: '8px',
            height: '8px',
            background: 'white',
            borderRadius: '50%',
            animation: 'pulse 1s infinite'
          }} />
          Recording...
        </div>
      )}

      <style>{`
        @keyframes pulse {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.5; }
        }
      `}</style>
    </div>
  );
}

// Export the recording functions for external use
export const useVideoRecorder = (targetRef: React.RefObject<HTMLElement>) => {
  const [isRecording, setIsRecording] = useState(false);
  const [recordingBlob, setRecordingBlob] = useState<Blob | null>(null);

  const startRecording = useCallback(async () => {
    // Implementation would be similar to above
    setIsRecording(true);
  }, []);

  const stopRecording = useCallback(async () => {
    setIsRecording(false);
    // Implementation would return blob
  }, []);

  return {
    isRecording,
    recordingBlob,
    startRecording,
    stopRecording
  };
}; 