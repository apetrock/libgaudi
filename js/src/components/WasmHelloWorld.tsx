import { useState, useRef } from 'react';
import { useWasmModule, WasmModule } from '../utils/wasmLoader';

/**
 * Phase 1: Basic WASM "Hello World" Integration Component
 * 
 * This component demonstrates loading and calling C++ functions from WebAssembly
 * before moving on to more complex geometry logging integration.
 */

interface HelloModule extends WasmModule {
  hello_world(): string;
  add_numbers(a: number, b: number): number;
  multiply_floats(x: number, y: number): number;
  // Logger functions
  log_line(x0: number, y0: number, z0: number, x1: number, y1: number, z1: number, 
           r: number, g: number, b: number, a: number): void;
  clear_logger(): void;
  test_logger(): string;
  get_line_count(): number;
  get_logger_status(): string;
  MathHelper: new (initialValue: number) => {
    getValue(): number;
    setValue(value: number): void;
    calculate(input: number): number;
    getStatus(): string;
  };
}

export function WasmHelloWorld() {
  const { loading, error, module } = useWasmModule<HelloModule, void>('hello_world');
  
  const [results, setResults] = useState({
    helloMessage: '',
    addResult: 0,
    multiplyResult: 0,
    mathHelperStatus: '',
    mathHelperCalculation: 0,
    // Logger results
    loggerTestResult: '',
    loggerStatus: '',
    lineCount: 0
  });
  
  const [inputs, setInputs] = useState({
    addA: 5,
    addB: 7,
    multiplyX: 3.14,
    multiplyY: 2.0,
    mathHelperValue: 10.0,
    mathHelperInput: 5.0
  });
  
  const mathHelperRef = useRef<any>(null);

  const runTests = (module: HelloModule) => {
    try {
      const helloMessage = module.hello_world();
      const addResult = module.add_numbers(inputs.addA, inputs.addB);
      const multiplyResult = module.multiply_floats(inputs.multiplyX, inputs.multiplyY);
      
      let mathHelperStatus = '';
      let mathHelperCalculation = 0;
      
      if (mathHelperRef.current) {
        mathHelperStatus = mathHelperRef.current.getStatus();
        mathHelperCalculation = mathHelperRef.current.calculate(inputs.mathHelperInput);
      }
      
      // Test logger functions
      const loggerTestResult = module.test_logger();
      
      // Log a few more lines to test the logger
      module.log_line(0, 0, 0, 2, 0, 0, 0, 1, 0, 1); // Green line along X axis
      module.log_line(0, 0, 0, 0, 2, 0, 0, 0, 1, 1); // Blue line along Y axis
      module.log_line(0, 0, 0, 0, 0, 2, 1, 1, 0, 1); // Yellow line along Z axis
      
      const loggerStatus = module.get_logger_status();
      const lineCount = module.get_line_count();
      
      setResults({
        helloMessage,
        addResult,
        multiplyResult,
        mathHelperStatus,
        mathHelperCalculation,
        loggerTestResult,
        loggerStatus,
        lineCount
      });
    } catch (error) {
      console.error('Error running WASM tests:', error);
    }
  };

  // Initialize MathHelper when module loads
  if (module && !mathHelperRef.current) {
    mathHelperRef.current = new module.MathHelper(inputs.mathHelperValue);
    // Run initial tests
    runTests(module);
  }

  const handleRunTests = () => {
    if (module && mathHelperRef.current) {
      // Update MathHelper value if changed
      mathHelperRef.current.setValue(inputs.mathHelperValue);
      runTests(module);
    }
  };

  const handleClearLogger = () => {
    if (module) {
      module.clear_logger();
      const loggerStatus = module.get_logger_status();
      const lineCount = module.get_line_count();
      
      setResults(prev => ({
        ...prev,
        loggerStatus,
        lineCount
      }));
    }
  };

  if (loading) {
    return (
      <div style={{ 
        width: '100%', 
        height: '100%', 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        background: '#1a1a2e',
        color: 'white',
        fontFamily: 'monospace'
      }}>
        <div style={{ textAlign: 'center' }}>
          <div style={{ 
            width: '40px', 
            height: '40px', 
            border: '4px solid rgba(74, 144, 226, 0.3)',
            borderTop: '4px solid #4a90e2',
            borderRadius: '50%',
            animation: 'spin 1s linear infinite',
            margin: '0 auto 20px'
          }} />
          <h3 style={{ color: '#4a90e2' }}>Loading WebAssembly Module...</h3>
          <p style={{ color: '#ccc' }}>Initializing C++ Hello World</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div style={{ 
        width: '100%', 
        height: '100%', 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        background: '#1a1a2e',
        color: 'white',
        fontFamily: 'monospace'
      }}>
        <div style={{ textAlign: 'center', maxWidth: '600px' }}>
          <h3 style={{ color: '#e74c3c', marginBottom: '20px' }}>❌ WASM Loading Failed</h3>
          <div style={{ 
            background: 'rgba(231, 76, 60, 0.1)',
            border: '1px solid rgba(231, 76, 60, 0.3)',
            borderRadius: '8px',
            padding: '20px',
            marginBottom: '20px'
          }}>
            <p style={{ margin: '0', color: '#ccc' }}>{error}</p>
          </div>
          
          <div style={{ 
            background: 'rgba(74, 144, 226, 0.1)',
            border: '1px solid rgba(74, 144, 226, 0.3)',
            borderRadius: '8px',
            padding: '20px',
            textAlign: 'left'
          }}>
            <h4 style={{ color: '#4a90e2', margin: '0 0 15px 0' }}>To build WASM module:</h4>
            <ol style={{ color: '#ccc', lineHeight: '1.6' }}>
              <li>Install Emscripten SDK</li>
              <li>Navigate to: <code>wasm/</code> directory</li>
              <li>Run: <code>./build.ps1</code> (Windows) or <code>make</code> (Linux/Mac)</li>
              <li>Refresh this page</li>
            </ol>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div style={{ 
      width: '100%', 
      height: '100%', 
      padding: '80px 20px 20px',
      background: '#1a1a2e',
      color: 'white',
      fontFamily: 'monospace',
      overflow: 'auto'
    }}>
      <div style={{ maxWidth: '800px', margin: '0 auto' }}>
        <h2 style={{ color: '#4a90e2', textAlign: 'center', marginBottom: '30px' }}>
          🚀 WASM Hello World - Phase 1
        </h2>
        
        {/* Status */}
        <div style={{ 
          background: 'rgba(46, 204, 113, 0.1)',
          border: '1px solid rgba(46, 204, 113, 0.3)',
          borderRadius: '8px',
          padding: '15px',
          marginBottom: '30px',
          textAlign: 'center'
        }}>
          <h3 style={{ color: '#2ecc71', margin: '0 0 10px 0' }}>✅ WebAssembly Module Loaded</h3>
          <p style={{ margin: '0', color: '#ccc' }}>C++ functions are ready to call from JavaScript</p>
        </div>

        {/* Results */}
        <div style={{ 
          background: 'rgba(74, 144, 226, 0.1)',
          border: '1px solid rgba(74, 144, 226, 0.3)',
          borderRadius: '8px',
          padding: '20px',
          marginBottom: '30px'
        }}>
          <h3 style={{ color: '#4a90e2', margin: '0 0 15px 0' }}>🎯 Function Results</h3>
          
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '15px' }}>
            <div>
              <strong style={{ color: '#2ecc71' }}>Hello World:</strong>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
                {results.helloMessage}
              </div>
            </div>
            
            <div>
              <strong style={{ color: '#2ecc71' }}>Add Numbers:</strong>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
                {inputs.addA} + {inputs.addB} = {results.addResult}
              </div>
            </div>
            
            <div>
              <strong style={{ color: '#2ecc71' }}>Multiply Floats:</strong>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
                {inputs.multiplyX} × {inputs.multiplyY} = {results.multiplyResult.toFixed(3)}
              </div>
            </div>
            
            <div>
              <strong style={{ color: '#2ecc71' }}>Math Helper:</strong>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
                Result: {results.mathHelperCalculation.toFixed(2)}
              </div>
            </div>
          </div>
            <div style={{ marginTop: '15px' }}>
            <strong style={{ color: '#2ecc71' }}>Math Helper Status:</strong>
            <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
              {results.mathHelperStatus}
            </div>
          </div>

          {/* Logger Results */}
          <div style={{ marginTop: '15px' }}>
            <strong style={{ color: '#e67e22' }}>📊 Line Logger:</strong>
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '10px', marginTop: '10px' }}>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px' }}>
                <strong>Lines Count:</strong> {results.lineCount}
              </div>
              <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px' }}>
                <strong>Status:</strong> {results.loggerStatus}
              </div>
            </div>
            <div style={{ background: 'rgba(0,0,0,0.3)', padding: '8px', borderRadius: '4px', marginTop: '5px' }}>
              <strong>Test Result:</strong> {results.loggerTestResult}
            </div>
          </div>
        </div>

        {/* Interactive Controls */}
        <div style={{ 
          background: 'rgba(74, 144, 226, 0.1)',
          border: '1px solid rgba(74, 144, 226, 0.3)',
          borderRadius: '8px',
          padding: '20px'
        }}>
          <h3 style={{ color: '#4a90e2', margin: '0 0 15px 0' }}>🎮 Interactive Controls</h3>
          
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '20px', marginBottom: '20px' }}>
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Add A:</label>
              <input 
                type="number" 
                value={inputs.addA}
                onChange={(e) => setInputs({...inputs, addA: parseInt(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
            
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Add B:</label>
              <input 
                type="number" 
                value={inputs.addB}
                onChange={(e) => setInputs({...inputs, addB: parseInt(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
            
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Multiply X:</label>
              <input 
                type="number" 
                step="0.01"
                value={inputs.multiplyX}
                onChange={(e) => setInputs({...inputs, multiplyX: parseFloat(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
            
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Multiply Y:</label>
              <input 
                type="number" 
                step="0.01"
                value={inputs.multiplyY}
                onChange={(e) => setInputs({...inputs, multiplyY: parseFloat(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
            
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Math Helper Value:</label>
              <input 
                type="number" 
                step="0.1"
                value={inputs.mathHelperValue}
                onChange={(e) => setInputs({...inputs, mathHelperValue: parseFloat(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
            
            <div>
              <label style={{ display: 'block', marginBottom: '5px', color: '#ccc' }}>Math Helper Input:</label>
              <input 
                type="number" 
                step="0.1"
                value={inputs.mathHelperInput}
                onChange={(e) => setInputs({...inputs, mathHelperInput: parseFloat(e.target.value) || 0})}
                style={{ width: '100%', padding: '8px', borderRadius: '4px', border: 'none', background: 'rgba(0,0,0,0.3)', color: 'white' }}
              />
            </div>
          </div>
            <button 
            onClick={handleRunTests}
            style={{
              width: '100%',
              padding: '12px',
              background: '#4a90e2',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontSize: '16px',
              fontWeight: 'bold',
              marginBottom: '10px'
            }}
          >
            🔄 Run Tests
          </button>
          
          <button 
            onClick={handleClearLogger}
            style={{
              width: '100%',
              padding: '12px',
              background: '#e67e22',
              color: 'white',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontSize: '16px',
              fontWeight: 'bold'
            }}
          >
            🗑️ Clear Logger
          </button>
        </div>
      </div>
    </div>
  );
}
