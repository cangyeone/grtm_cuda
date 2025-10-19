cmake_minimum_required(VERSION 3.12)
project(GRTM_CUDA LANGUAGES CXX CUDA)
set(CMAKE_CUDA_ARCHITECTURES OFF)
# Sets the possible install locations
set(CUBLAS_HINTS
  ${CUDA_ROOT}
  $ENV{CUDA_ROOT}
  $ENV{CUDA_TOOLKIT_ROOT_DIR}
)
set(CUBLAS_PATHS
  /usr
  /usr/local
  /usr/local/cuda
)



# Finds the include directories
find_path(CUBLAS_INCLUDE_DIRS
  NAMES cublas_v2.h cuda.h
  HINTS ${CUBLAS_HINTS}
  PATH_SUFFIXES include inc include/x86_64 include/x64
  PATHS ${CUBLAS_PATHS}
  DOC "cuBLAS include header cublas_v2.h"
)
mark_as_advanced(CUBLAS_INCLUDE_DIRS)

# Finds the libraries
find_library(CUDA_LIBRARIES
  NAMES cudart
  HINTS ${CUBLAS_HINTS}
  PATH_SUFFIXES lib lib64 lib/x86_64 lib/x64 lib/x86 lib/Win32 lib/import lib64/import
  PATHS ${CUBLAS_PATHS}
  DOC "CUDA library"
)
mark_as_advanced(CUDA_LIBRARIES)
find_library(CUBLAS_LIBRARIES
  NAMES cublas
  HINTS ${CUBLAS_HINTS}
  PATH_SUFFIXES lib lib64 lib/x86_64 lib/x64 lib/x86 lib/Win32 lib/import lib64/import
  PATHS ${CUBLAS_PATHS}
  DOC "cuBLAS library"
)
mark_as_advanced(CUBLAS_LIBRARIES)

find_library(CUFFT_LIBRARIES
  NAMES cufft
  HINTS ${CUBLAS_HINTS}
  PATH_SUFFIXES lib lib64 lib/x86_64 lib/x64 lib/x86 lib/Win32 lib/import lib64/import
  PATHS ${CUBLAS_PATHS}
  DOC "cufft library"
)
mark_as_advanced(CUFFT_LIBRARIES)

# ==================================================================================================

# Notification messages
if(NOT CUBLAS_INCLUDE_DIRS)
    message(STATUS "Could NOT find 'cuBLAS.h', install CUDA/cuBLAS or set CUDA_ROOT")
endif()
if(NOT CUDA_LIBRARIES)
    message(STATUS "Could NOT find CUDA library, install it or set CUDA_ROOT")
endif()
if(NOT CUBLAS_LIBRARIES)
    message(STATUS "Could NOT find cuBLAS library, install it or set CUDA_ROOT")
endif()

# Determines whether or not cuBLAS was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cuBLAS DEFAULT_MSG CUBLAS_INCLUDE_DIRS CUDA_LIBRARIES CUBLAS_LIBRARIES CUFFT_LIBRARIES)

# ==================================================================================================
file(GLOB_RECURSE POPM_SRC cusrc/popm/func/*.cu)
add_library(popm SHARED ${POPM_SRC})

file(GLOB_RECURSE SCILIB_SRC cusrc/scilib/func/*.cu)
add_library(scilib SHARED ${SCILIB_SRC})
target_link_libraries(popm scilib cufft cublas cuda)

#set_property(TARGET popm PROPERTY CUDA_ARCHITECTURES 70)
#set_property(TARGET scilib PROPERTY CUDA_ARCHITECTURES 70)

set(SOURCES cusrc/main.cu)
add_executable(test ${SOURCES})
target_compile_features(test PRIVATE cxx_std_11)
target_link_libraries(test PRIVATE ${CUBLAS_LIBRARIES} popm)
#set_property(TARGET test PROPERTY CUDA_ARCHITECTURES 70)