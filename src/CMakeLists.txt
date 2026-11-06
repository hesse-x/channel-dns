if(CHANNEL_DNS_FFT_BACKEND STREQUAL "FFTW")
 set(fft_backend fft_backend_fftw.f90)
 set(fft_libraries PkgConfig::FFTW)
else()
 set(fft_backend fft_backend_fftpack.f90)
 add_library(dns_fftpack STATIC ${CMAKE_SOURCE_DIR}/legacy/src/fftpack.f90)
 if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(dns_fftpack PRIVATE -std=legacy -w)
 endif()
 set(fft_libraries dns_fftpack)
endif()
add_library(channel_dns_fortran STATIC
  kinds.f90 toml.f90 config.f90 mesh.f90 ${fft_backend} fft.f90 equation_state.f90
  dns_vector.f90 dns_equations.f90 dns_timestep.f90 initial.f90 io.f90 statistics.f90
)
set_target_properties(channel_dns_fortran PROPERTIES
  Fortran_STANDARD 2008 Fortran_STANDARD_REQUIRED ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(channel_dns_fortran PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_link_libraries(channel_dns_fortran PUBLIC MPI::MPI_Fortran ${fft_libraries})
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Component views must remain descriptors, never hidden full-field copies.
  set_source_files_properties(dns_vector.f90 dns_timestep.f90 PROPERTIES
    COMPILE_OPTIONS "-Warray-temporaries;-Werror=array-temporaries")
  target_compile_options(channel_dns_fortran PRIVATE -Wall -Wextra -Wno-compare-reals -Wno-aliasing -ffree-line-length-none
    $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> $<$<CONFIG:Debug>:-fbacktrace>)
endif()

add_executable(channel_dns main.f90)
target_link_libraries(channel_dns PRIVATE channel_dns_fortran)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
 target_compile_options(channel_dns PRIVATE -ffree-line-length-none $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> -fbacktrace)
endif()

if(ENABLE_UNIT_TESTS)
  find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter)
  add_executable(fortran_numerics ${CMAKE_SOURCE_DIR}/tests/fortran/test_numerics.f90)
  target_link_libraries(fortran_numerics PRIVATE channel_dns_fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fortran_numerics PRIVATE -ffree-line-length-none
      $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> -fbacktrace)
  endif()
  foreach(ranks IN ITEMS 1 2 3 4)
    set(run_dir ${CMAKE_CURRENT_BINARY_DIR}/test-np${ranks})
    file(MAKE_DIRECTORY ${run_dir})
    add_test(NAME fortran_numerics_np${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:fortran_numerics> ${MPIEXEC_POSTFLAGS} ${run_dir})
  endforeach()
  add_test(NAME fortran_io_serial COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/tools/test_io_compatibility.py $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR})
  add_test(NAME fortran_io_mpi COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/tools/test_io_compatibility.py $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR}
    --mpiexec ${MPIEXEC_EXECUTABLE})
  add_test(NAME fortran_migration COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/tests/tools/test_fortran_migration.py $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR}
    --mpiexec ${MPIEXEC_EXECUTABLE})
  set_tests_properties(fortran_numerics_np1 fortran_numerics_np2 fortran_numerics_np3 fortran_numerics_np4
    fortran_io_serial fortran_io_mpi fortran_migration PROPERTIES
    TIMEOUT 180 ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endif()

if(ENABLE_UNIT_TESTS)
  add_executable(fortran_config ${CMAKE_SOURCE_DIR}/tests/fortran/test_config.f90)
  add_executable(fortran_zero_mode ${CMAKE_SOURCE_DIR}/tests/fortran/test_zero_mode.f90)
  foreach(target IN ITEMS fortran_config fortran_zero_mode)
    target_link_libraries(${target} PRIVATE channel_dns_fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${target} PRIVATE -ffree-line-length-none
        $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> -fbacktrace)
    endif()
  endforeach()
  add_test(NAME fortran_config COMMAND fortran_config ${CMAKE_SOURCE_DIR}/tests/fixtures/toml ${CMAKE_SOURCE_DIR})
  set_tests_properties(fortran_config PROPERTIES TIMEOUT 30)
  foreach(ranks IN ITEMS 1 2 3 4)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zero-mode-np${ranks})
    add_test(NAME fortran_zero_mode_np${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:fortran_zero_mode> ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/zero-mode-np${ranks})
    set_tests_properties(fortran_zero_mode_np${ranks} PROPERTIES TIMEOUT 60
      ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
  endforeach()
endif()

if(ENABLE_UNIT_TESTS)
  foreach(suite IN ITEMS field_io vector_mesh)
    add_executable(fortran_${suite} ${CMAKE_SOURCE_DIR}/tests/fortran/test_${suite}.f90)
    target_link_libraries(fortran_${suite} PRIVATE channel_dns_fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(fortran_${suite} PRIVATE -ffree-line-length-none
        $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> -fbacktrace)
    endif()
    foreach(ranks IN ITEMS 1 2 3 4)
      set(run_dir ${CMAKE_CURRENT_BINARY_DIR}/${suite}-np${ranks})
      file(MAKE_DIRECTORY ${run_dir})
      add_test(NAME fortran_${suite}_np${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:fortran_${suite}> ${MPIEXEC_POSTFLAGS} ${run_dir})
      set_tests_properties(fortran_${suite}_np${ranks} PROPERTIES TIMEOUT 60
        ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
    endforeach()
  endforeach()
endif()

if(ENABLE_UNIT_TESTS)
  add_executable(fortran_fused_solver ${CMAKE_SOURCE_DIR}/tests/fortran/test_fused_solver.f90)
  target_link_libraries(fortran_fused_solver PRIVATE channel_dns_fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fortran_fused_solver PRIVATE -ffree-line-length-none
      $<$<CONFIG:Debug>:-fcheck=bounds,do,mem,pointer> -fbacktrace)
  endif()
  foreach(ranks IN ITEMS 1 3 4)
    set(run_dir ${CMAKE_CURRENT_BINARY_DIR}/fused-np${ranks})
    file(MAKE_DIRECTORY ${run_dir})
    add_test(NAME fortran_fused_solver_np${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:fortran_fused_solver> ${MPIEXEC_POSTFLAGS} ${run_dir})
    set_tests_properties(fortran_fused_solver_np${ranks} PROPERTIES TIMEOUT 60
      ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
  endforeach()
endif()

if(CHANNEL_DNS_FFT_BACKEND STREQUAL "FFTW")
 target_include_directories(channel_dns_fortran PRIVATE ${FFTW_FORTRAN_INCLUDE_DIR})
endif()
