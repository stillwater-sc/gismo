### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## CoDiPack extension
project(gsCoDiPackExtension)

if(CMAKE_CXX_STANDARD EQUAL "98" )
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for CoDiPack extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

gismo_fetch_directory(CoDiPack
    URL https://github.com/SciCompKL/CoDiPack/archive/numericExtension.zip
    DESTINATION    external
)

set(CODIPACK_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/CoDiPack/include CACHE INTERNAL "")
include_directories(${CODIPACK_INCLUDE_DIR})

install(DIRECTORY ${CODIPACK_INCLUDE_DIR}
        DESTINATION include/codi/ OPTIONAL)
