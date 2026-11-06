find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter)
set(legacy_main ${CMAKE_CURRENT_SOURCE_DIR}/src/channel_with_Microblow_ycut.f90)
set(legacy_kernels
  src/EQSOLVE_channel_ycut.f90 src/fft_channel_ycut.f90 src/fftpack.f90
  src/for_channel_Microblow_ycut_.f90 src/for_channel_statistics_ycut.f90
  src/STEP23_channel_ycut.f90 src/vector_channel_ycut.f90)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/benchmark_main.f90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/instrument.py
    ${legacy_main} ${CMAKE_CURRENT_BINARY_DIR}/benchmark_main.f90
  DEPENDS ${legacy_main} ${CMAKE_CURRENT_SOURCE_DIR}/instrument.py VERBATIM)
add_executable(channel_dns_legacy ${CMAKE_CURRENT_BINARY_DIR}/benchmark_main.f90 ${legacy_kernels})
add_executable(channel_dns_legacy_raw EXCLUDE_FROM_ALL ${legacy_main} ${legacy_kernels})
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/aligned_equations.f90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/align_zero_mode.py
    ${CMAKE_CURRENT_SOURCE_DIR}/src/EQSOLVE_channel_ycut.f90 ${CMAKE_CURRENT_BINARY_DIR}/aligned_equations.f90
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/EQSOLVE_channel_ycut.f90 ${CMAKE_CURRENT_SOURCE_DIR}/align_zero_mode.py VERBATIM)
set(aligned_kernels ${legacy_kernels})
list(REMOVE_ITEM aligned_kernels src/EQSOLVE_channel_ycut.f90)
add_executable(channel_dns_legacy_aligned ${CMAKE_CURRENT_BINARY_DIR}/benchmark_main.f90
  ${CMAKE_CURRENT_BINARY_DIR}/aligned_equations.f90 ${aligned_kernels})
add_custom_target(legacy_generated DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/benchmark_main.f90
  ${CMAKE_CURRENT_BINARY_DIR}/aligned_equations.f90)
foreach(target IN ITEMS channel_dns_legacy channel_dns_legacy_raw channel_dns_legacy_aligned)
  add_dependencies(${target} legacy_generated)
  set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${target}-modules)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/${target}-modules)
  target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # FFTPACK uses legacy assumed-size conventions and MPI's mpif.h binding.
    target_compile_options(${target} PRIVATE -std=legacy -ffree-line-length-none -fallow-argument-mismatch -w)
  endif()
endforeach()

if(ENABLE_UNIT_TESTS)
  add_test(NAME legacy_precision COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/tools/test_legacy.py $<TARGET_FILE:channel_dns>
    $<TARGET_FILE:channel_dns_legacy> $<TARGET_FILE:channel_dns_legacy_aligned>
    ${MPIEXEC_EXECUTABLE} ${CMAKE_SOURCE_DIR})
  set_tests_properties(legacy_precision PROPERTIES TIMEOUT 180
    ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endif()
