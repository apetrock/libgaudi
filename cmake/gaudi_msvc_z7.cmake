# MSVC: Prefer /Z7 (debug info in .obj) over /Zi (shared vc*.pdb per project).
# C1041 happens when multiple cl.exe writers contend for vc143.pdb; /FS and
# UseMultiToolTask=false do not always fix it on VS 17+ with parallel builds.

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  return()
endif()

foreach(_lang C CXX)
  foreach(_suffix DEBUG RELWITHDEBINFO ASAN)
    set(_var CMAKE_${_lang}_FLAGS_${_suffix})
    if(DEFINED ${_var})
      string(REPLACE "/Zi" "/Z7" _fixed "${${_var}}")
      if(NOT _fixed STREQUAL "${${_var}}")
        set(${_var} "${_fixed}" CACHE STRING "Patched: /Z7 avoids MSVC C1041 on vc*.pdb" FORCE)
      endif()
    endif()
  endforeach()
endforeach()
