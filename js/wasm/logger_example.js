// Example of how to use the callback-based logging system

// Set up custom logging callbacks
Module.set_info_callback(function(message) {
    console.log('[INFO]', message);
    // You could also send this to a React component or other UI
    // window.dispatchEvent(new CustomEvent('log-info', { detail: message }));
});

Module.set_warning_callback(function(message) {
    console.warn('[WARNING]', message);
    // window.dispatchEvent(new CustomEvent('log-warning', { detail: message }));
});

Module.set_error_callback(function(message) {
    console.error('[ERROR]', message);
    // window.dispatchEvent(new CustomEvent('log-error', { detail: message }));
});

Module.set_debug_callback(function(message) {
    console.debug('[DEBUG]', message);
    // window.dispatchEvent(new CustomEvent('log-debug', { detail: message }));
});

// Example React component integration
/*
function LogDisplay() {
    const [logs, setLogs] = useState([]);
    
    useEffect(() => {
        const handleLog = (event) => {
            setLogs(prev => [...prev, { 
                type: event.type.replace('log-', ''), 
                message: event.detail,
                timestamp: new Date()
            }]);
        };
        
        window.addEventListener('log-info', handleLog);
        window.addEventListener('log-warning', handleLog);
        window.addEventListener('log-error', handleLog);
        window.addEventListener('log-debug', handleLog);
        
        return () => {
            window.removeEventListener('log-info', handleLog);
            window.removeEventListener('log-warning', handleLog);
            window.removeEventListener('log-error', handleLog);
            window.removeEventListener('log-debug', handleLog);
        };
    }, []);
    
    return (
        <div className="log-display">
            {logs.map((log, index) => (
                <div key={index} className={`log-entry log-${log.type}`}>
                    <span className="timestamp">{log.timestamp.toLocaleTimeString()}</span>
                    <span className="type">[{log.type.toUpperCase()}]</span>
                    <span className="message">{log.message}</span>
                </div>
            ))}
        </div>
    );
}
*/

// Now when C++ code calls gaudi::logger::info << "Hello World" << std::endl;
// it will be routed through the JavaScript callback and can be displayed in the UI

// Note: The callbacks are set up in the WASM terminal logger module,
// so they will be available as Module.set_info_callback, etc. 