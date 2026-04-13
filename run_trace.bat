@echo off
setlocal
set ROD=%~dp0build_gl_viewers\projects\rod_constraints_test\Release\rod_constraints_test.exe
set GROW=%~dp0build_gl_viewers\projects\growth_study\Release\growth_study.exe
set OUT=%~dp0_trace_capture.txt
echo === rod_constraints_test ===> "%OUT%"
"%ROD%" >> "%OUT%" 2>&1
echo EXIT_ROD=%ERRORLEVEL%>> "%OUT%"
echo === growth_study ===>> "%OUT%"
"%GROW%" >> "%OUT%" 2>&1
echo EXIT_GROW=%ERRORLEVEL%>> "%OUT%"
