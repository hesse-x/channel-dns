cmake_minimum_required(VERSION 3.16)
project(channel_dns CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CHANNEL_DNS_PURE_MPI "Execute CpuDevice inline on each rank; omit compute workers" OFF)

option(ENABLE_UNIT_TESTS "Build unit tests" ON)

find_package(MPI REQUIRED CXX)      # D7：MPI 必选
if(NOT CHANNEL_DNS_PURE_MPI)
  find_package(Threads REQUIRED)
  find_library(FFTW_THREADS_LIBRARY NAMES fftw3_threads)
  if(NOT FFTW_THREADS_LIBRARY)
    message(FATAL_ERROR "FFTW POSIX threads library fftw3_threads is required")
  endif()
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

add_compile_options(-Wall)
add_compile_options($<$<CONFIG:Release>:-O3>)

add_subdirectory(src)
if(ENABLE_UNIT_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
