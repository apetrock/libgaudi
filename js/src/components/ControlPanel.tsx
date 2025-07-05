import React from 'react';
import type { ControlPanelProps } from '../types/simulation';

/**
 * Control panel component for adjusting simulation parameters
 */
export function ControlPanel({
  params,
  onParamsChange,
  isPlaying,
  onPlayPause,
  onReset,
  onStep
}: ControlPanelProps) {
  const handleParamChange = (key: keyof typeof params, value: number) => {
    onParamsChange({ [key]: value });
  };

  return (
    <div className="rod-simulation-controls" style={{
      position: 'absolute',
      top: 10,
      left: 10,
      background: 'rgba(0, 0, 0, 0.8)',
      color: 'white',
      padding: '15px',
      borderRadius: '8px',
      fontFamily: 'monospace',
      fontSize: '12px',
      zIndex: 1000,
      minWidth: '250px'
    }}>
      <h3 style={{ margin: '0 0 15px 0' }}>Rod Simulation Controls</h3>
      
      {/* Playback Controls */}
      <div style={{ marginBottom: '15px' }}>
        <button 
          onClick={onPlayPause}
          style={{
            marginRight: '8px',
            padding: '5px 10px',
            background: isPlaying ? '#e74c3c' : '#2ecc71',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer'
          }}
        >
          {isPlaying ? 'Pause' : 'Play'}
        </button>
        
        <button 
          onClick={onStep}
          style={{
            marginRight: '8px',
            padding: '5px 10px',
            background: '#3498db',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer'
          }}
        >
          Step
        </button>
        
        <button 
          onClick={onReset}
          style={{
            padding: '5px 10px',
            background: '#f39c12',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer'
          }}
        >
          Reset
        </button>
      </div>

      {/* Parameter Controls */}
      <div>
        <div style={{ marginBottom: '10px' }}>
          <label>Growth Rate: {params.growthRate.toFixed(2)}</label>
          <input
            type="range"
            min="0.1"
            max="3.0"
            step="0.1"
            value={params.growthRate}
            onChange={(e) => handleParamChange('growthRate', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>Constraint Strength: {params.constraintStrength.toFixed(3)}</label>
          <input
            type="range"
            min="0.001"
            max="0.5"
            step="0.001"
            value={params.constraintStrength}
            onChange={(e) => handleParamChange('constraintStrength', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>Collision Strength: {params.collisionStrength.toFixed(1)}</label>
          <input
            type="range"
            min="0.1"
            max="5.0"
            step="0.1"
            value={params.collisionStrength}
            onChange={(e) => handleParamChange('collisionStrength', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>Time Step: {params.timeStep.toFixed(3)}</label>
          <input
            type="range"
            min="0.001"
            max="0.1"
            step="0.001"
            value={params.timeStep}
            onChange={(e) => handleParamChange('timeStep', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>Bending Stiffness: {params.bendingStiffness.toFixed(3)}</label>
          <input
            type="range"
            min="0.001"
            max="0.2"
            step="0.001"
            value={params.bendingStiffness}
            onChange={(e) => handleParamChange('bendingStiffness', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>Stretch Stiffness: {params.stretchStiffness.toFixed(3)}</label>
          <input
            type="range"
            min="0.001"
            max="0.2"
            step="0.001"
            value={params.stretchStiffness}
            onChange={(e) => handleParamChange('stretchStiffness', parseFloat(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>

        <div style={{ marginBottom: '10px' }}>
          <label>SDF Transition Frame: {params.sdfTransitionFrame}</label>
          <input
            type="range"
            min="100"
            max="1000"
            step="50"
            value={params.sdfTransitionFrame}
            onChange={(e) => handleParamChange('sdfTransitionFrame', parseInt(e.target.value))}
            style={{ width: '100%', marginTop: '2px' }}
          />
        </div>
      </div>
    </div>
  );
}
