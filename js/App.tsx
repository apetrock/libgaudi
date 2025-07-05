import React from 'react';
import { BrowserRouter as Router, Routes, Route, Link, useLocation } from 'react-router-dom';
import { BoxDemo } from './src/components/BoxDemo';
import { ZustandLineLoggerTest } from './src/components/ZustandLineLoggerTest';
import { WasmHelloWorld } from './src/components/WasmHelloWorld';
import { GaudiLoggerTest } from './src/components/GaudiLoggerTest';
import { RodConstraintsTest } from './src/components/RodConstraintsTest';
import { Button } from './src/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './src/components/ui/card';
import { Badge } from './src/components/ui/badge';
import { 
  Box, 
  TestTube, 
  Waves, 
  Microscope, 
  Rocket, 
  ArrowRight,
  Github,
  Play,
  Pause
} from 'lucide-react';
import { Foo_demo } from './src/components/Foo_demo';
import { Bar_demo } from './src/components/Bar_demo';
/**
 * Navigation component for different test scenarios
 */
function Navigation() {
  const location = useLocation();
    const navItems = [
    { path: '/box-demo', label: 'Box Demo', icon: Box, status: 'stable' },
    { path: '/line-logger', label: 'Line Logger', icon: Waves, status: 'stable' },
    { path: '/gaudi-logger', label: 'Gaudi Logger', icon: Microscope, status: 'beta' },
    { path: '/rod-constraints', label: 'Rod Constraints', icon: TestTube, status: 'beta' },
    { path: '/wasm-hello', label: 'WASM Hello', icon: TestTube, status: 'stable' },
    { path: '/future-test', label: 'Future', icon: Rocket, status: 'planned' },
  ];

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'stable': return 'default';
      case 'beta': return 'secondary';
      case 'planned': return 'outline';
      default: return 'default';
    }
  };

  return (
    <nav className="fixed top-0 left-0 right-0 z-50 bg-background/80 backdrop-blur-sm border-b border-gray-200">
      <div className="container mx-auto px-4 py-3">
        <div className="flex items-center justify-between">          <div className="flex items-center space-x-4">
            <Link to="/" className="flex items-center space-x-2 hover:opacity-80 transition-opacity">
              <Microscope className="h-6 w-6 text-primary" />
              <h1 className="text-lg font-bold text-foreground">
                Rod Simulation
              </h1>
              <Badge variant="outline" className="text-xs">
                Phase 1.5
              </Badge>
            </Link>
          </div>          <div className="flex flex-row items-center space-x-2">
            {navItems.map((item) => {
              const Icon = item.icon;
              const isActive = location.pathname === item.path;
              
              return (
                <Button
                  key={item.path}
                  variant={isActive ? "default" : "ghost"}
                  size="sm"
                  asChild
                  className="relative"
                >
                  <Link to={item.path} className="flex items-center space-x-2">
                    <Icon className="h-4 w-4" />
                    <span className="hidden md:inline">{item.label}</span>
                  </Link>
                </Button>
              );
            })}
          </div>
          
          <div className="flex items-center space-x-2">
            <Button variant="ghost" size="icon" asChild>
              <a 
                href="https://github.com/your-repo/rod-simulation" 
                target="_blank" 
                rel="noopener noreferrer"
                className="text-muted-foreground hover:text-foreground"
              >
                <Github className="h-4 w-4" />
              </a>
            </Button>
          </div>
        </div>
      </div>
    </nav>
  );
}

/**
 * Box Demo page - Phase 0 milestone validation
 */
function BoxDemoPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4" style={{ paddingTop: '5rem' }}>
        <div className="h-[calc(100vh-5rem)]">
          <BoxDemo 
            color="#4a90e2"
            rotationSpeed={1.0}
            showControls={true}
          />
        </div>
      </div>
    </div>
  );
}

/**
 * Line Logger Test page - validates debug line rendering
 */
function LineLoggerPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4" style={{ paddingTop: '5rem' }}>
        <div className="h-[calc(100vh-5rem)]">
          <ZustandLineLoggerTest />
        </div>
      </div>
    </div>
  );
}

/**
 * WASM Hello World page - Phase 1 basic integration
 */
function WasmHelloPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4" style={{ paddingTop: '5rem' }}>
        <div className="h-[calc(100vh-5rem)]">
          <WasmHelloWorld />
        </div>
      </div>
    </div>
  );
}

/**
 * Gaudi Logger Test page - Phase 1.5 gaudi::logger interface testing
 */
function GaudiLoggerPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4" style={{ paddingTop: '5rem' }}>
        <div className="h-[calc(100vh-5rem)]">
          <GaudiLoggerTest />
        </div>
      </div>
    </div>
  );
}

/**
 * Rod Constraints Test page - Phase 1.5 rod dynamics simulation
 */
function RodConstraintsPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4" style={{ paddingTop: '5rem' }}>
        <div className="h-[calc(100vh-5rem)]">
          <RodConstraintsTest />
        </div>
      </div>
    </div>
  );
}

/**
 * Future integration test page with beautiful design
 */
function FutureTestPage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4 py-8" style={{ paddingTop: '5rem' }}>
        <div className="max-w-4xl mx-auto">
          <div className="text-center mb-8">
            <div className="flex items-center justify-center space-x-2 mb-4">
              <Rocket className="h-8 w-8 text-primary" />
              <h1 className="text-4xl font-bold text-foreground">
                Phase 2: Full Integration
              </h1>
            </div>
            <p className="text-xl text-muted-foreground">
              Advanced WASM integration with complete rod physics simulation
            </p>
          </div>

          <div className="grid md:grid-cols-2 gap-6 mb-8">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <TestTube className="h-5 w-5 text-primary" />
                  <span>Current Progress</span>
                </CardTitle>
                <CardDescription>
                  Phase 1.5 milestones completed
                </CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <span className="text-sm">✅ Basic Box Rendering</span>
                    <Badge variant="default">Done</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">✅ Line Logger System</span>
                    <Badge variant="default">Done</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">✅ WASM Integration</span>
                    <Badge variant="default">Done</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">🔬 Gaudi Logger Interface</span>
                    <Badge variant="secondary">Testing</Badge>
                  </div>
                </div>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <ArrowRight className="h-5 w-5 text-primary" />
                  <span>Next Steps</span>
                </CardTitle>
                <CardDescription>
                  Phase 2 roadmap and goals
                </CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <span className="text-sm">🎯 Rod Physics Engine</span>
                    <Badge variant="outline">Planned</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">🎨 Advanced Visualization</span>
                    <Badge variant="outline">Design</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">⚡ Performance Optimization</span>
                    <Badge variant="outline">Research</Badge>
                  </div>
                  <div className="flex items-center justify-between">
                    <span className="text-sm">🎮 Interactive Controls</span>
                    <Badge variant="outline">Planning</Badge>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>

          <Card className="bg-gradient-to-r from-primary/10 to-secondary/10 border-primary/20">
            <CardHeader>
              <CardTitle className="text-2xl text-center">
                🚀 Rod Physics Simulation
              </CardTitle>
              <CardDescription className="text-center text-lg">
                Interactive browser-based physics simulation with real-time visualization
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="grid md:grid-cols-3 gap-4 text-center">
                <div className="space-y-2">
                  <div className="text-3xl">🌊</div>
                  <h3 className="font-semibold">Dynamic Constraints</h3>
                  <p className="text-sm text-muted-foreground">
                    Stretch, bend, twist, and collision constraints with real-time solving
                  </p>
                </div>
                <div className="space-y-2">
                  <div className="text-3xl">🎯</div>
                  <h3 className="font-semibold">Growth Mechanics</h3>
                  <p className="text-sm text-muted-foreground">
                    Adaptive rod growth based on SDF interactions and environmental factors
                  </p>
                </div>
                <div className="space-y-2">
                  <div className="text-3xl">✨</div>
                  <h3 className="font-semibold">Visual Excellence</h3>
                  <p className="text-sm text-muted-foreground">
                    Smooth tube geometry, color-coded growth, and interactive 3D controls
                  </p>
                </div>
              </div>
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}

/**
 * Home page - landing with overview
 */
function HomePage() {
  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto max-w-4xl px-4 py-8" style={{ paddingTop: '5rem' }}>
        <div className="max-w-4xl mx-auto">
          <div className="text-center mb-12">
            <div className="flex items-center justify-center space-x-3 mb-6">
              <Microscope className="h-12 w-12 text-primary" />
              <h1 className="text-5xl font-bold text-foreground">
                Rod Simulation
              </h1>
            </div>
            <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
              Interactive browser-based physics simulation powered by WebAssembly, 
              React Three Fiber, and advanced constraint solving
            </p>
            <div className="flex items-center justify-center space-x-4 mt-6">
              <Button size="lg" asChild>
                <Link to="/box-demo" className="flex items-center space-x-2">
                  <Play className="h-4 w-4" />
                  <span>Start Demo</span>
                </Link>
              </Button>
              <Button variant="outline" size="lg" asChild>
                <Link to="/gaudi-logger" className="flex items-center space-x-2">
                  <Microscope className="h-4 w-4" />
                  <span>Try Logger Test</span>
                </Link>
              </Button>
            </div>
          </div>

          <div className="grid md:grid-cols-3 gap-6">
            <Card className="hover:shadow-lg transition-shadow">
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <Box className="h-5 w-5 text-primary" />
                  <span>3D Visualization</span>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Real-time 3D rendering with React Three Fiber, smooth animations, 
                  and interactive orbital controls for immersive physics exploration.
                </p>
              </CardContent>
            </Card>

            <Card className="hover:shadow-lg transition-shadow">
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <TestTube className="h-5 w-5 text-primary" />
                  <span>WASM Performance</span>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  High-performance physics computation using WebAssembly, 
                  enabling complex simulations to run smoothly in the browser.
                </p>
              </CardContent>
            </Card>

            <Card className="hover:shadow-lg transition-shadow">
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <Waves className="h-5 w-5 text-primary" />
                  <span>Debug Visualization</span>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground">
                  Advanced debug line rendering system with real-time data streaming 
                  from C++ physics engine to JavaScript visualization layer.
                </p>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </div>
  );
}

/**
 * Demo application with React Router navigation between test scenarios
 */
function App() {
  return (
    <Router>
      <div className="min-h-screen bg-background">
        <Navigation />
        <Routes>
          <Route path="/" element={<HomePage />} />
          <Route path="/box-demo" element={<BoxDemoPage />} />
          <Route path="/line-logger" element={<LineLoggerPage />} />
          <Route path="/gaudi-logger" element={<GaudiLoggerPage />} />
          <Route path="/wasm-hello" element={<WasmHelloPage />} />
          <Route path="/future-test" element={<FutureTestPage />} />
          <Route path="/rod-constraints" element={<RodConstraintsPage />} />
          <Route path="/foo_demo" element={<Foo_demo />} />
          <Route path="/bar_demo" element={<Bar_demo />} />
        </Routes>
      </div>
    </Router>
  );
}

export default App;
