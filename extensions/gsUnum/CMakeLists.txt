### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Unum extension
project(gsUnumExtension)

if(CMAKE_CXX_STANDARD EQUAL "98")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for Unum extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

gismo_fetch_directory(Unum
    GIT_REPOSITORY https://github.com/stillwater-sc/universal.git
    DESTINATION    external
)

set(UNUM_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/Unum/include CACHE INTERNAL "")
include_directories(${UNUM_INCLUDE_DIR})

install(DIRECTORY ${UNUM_INCLUDE_DIR}
        DESTINATION include/unum/ OPTIONAL)
