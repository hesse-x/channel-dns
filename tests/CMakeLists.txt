include(FetchContent)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1)
# Windows  CRT 设置对我们无用，但保留官方建议开关以免警告
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(UNIT_TEST_SRCS
  device_test_main.cpp
  test_cpu_device.cpp
  test_array4d.cpp
  test_mesh.cpp
  test_params.cpp
  test_toml.cpp
  test_interfaces.cpp
  test_fft_x.cpp
  test_fft_z.cpp
  test_fft_y.cpp
  test_xyzfft.cpp
  test_vector_ops.cpp
  test_nonlinear_ops.cpp
  test_eqsolve.cpp
  test_timestepper.cpp
  test_initial_bc.cpp
  test_field_io.cpp
  test_statistics.cpp
  test_diagnostics.cpp
  test_simulation.cpp
)
if(CHANNEL_DNS_MPI_ONLY)
  list(APPEND UNIT_TEST_SRCS test_device_inline.cpp)
else()
  list(APPEND UNIT_TEST_SRCS test_device_queue.cpp)
endif()
add_executable(unit_tests ${UNIT_TEST_SRCS})
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(unit_tests PRIVATE CHANNEL_DNS_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(unit_tests PRIVATE
  GTest::gtest
  MPI::MPI_CXX
  channel_dns_core
  channel_dns_fft
  channel_dns_solver
  channel_dns_timestep
  channel_dns_io
  channel_dns_sim
)

include(GoogleTest)
gtest_discover_tests(unit_tests)

add_executable(mpi_transform_tests test_mpi_transforms.cpp)
target_include_directories(mpi_transform_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(mpi_transform_tests PRIVATE MPI::MPI_CXX channel_dns_core channel_dns_fft)
foreach(_ranks IN ITEMS 2 3 4)
  add_test(NAME mpi_transforms_np${_ranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_transform_tests> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_transforms_np${_ranks} PROPERTIES
    TIMEOUT 60
    ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endforeach()

find_package(Python3 3.11 COMPONENTS Interpreter REQUIRED)
add_test(NAME io_compatibility
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_io_compatibility.py
    $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR})
add_test(NAME io_compatibility_mpi
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_io_compatibility.py
    $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR} --mpiexec ${MPIEXEC_EXECUTABLE})
set_tests_properties(io_compatibility io_compatibility_mpi PROPERTIES
  TIMEOUT 120 ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")

add_executable(mpi_context_tests test_mpi_context.cpp)
target_link_libraries(mpi_context_tests PRIVATE GTest::gtest channel_dns_core)
foreach(_ranks IN ITEMS 1 2 3 4)
  add_test(NAME mpi_context_np${_ranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_context_tests> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_context_np${_ranks} PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endforeach()

if(NOT CHANNEL_DNS_MPI_ONLY)
# Exercise the same numerical suite with one persistent pool of eight threads.
add_test(NAME hybrid_unit_tests COMMAND $<TARGET_FILE:unit_tests> --threads 8)
set_tests_properties(hybrid_unit_tests PROPERTIES TIMEOUT 120)
foreach(_ranks IN ITEMS 2 3 4)
  add_test(NAME hybrid_transforms_np${_ranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_transform_tests> 4 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(hybrid_transforms_np${_ranks} PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endforeach()
add_test(NAME hybrid_io_compatibility_mpi
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_io_compatibility.py
    $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR} --mpiexec ${MPIEXEC_EXECUTABLE} --threads 4)
set_tests_properties(hybrid_io_compatibility_mpi PROPERTIES TIMEOUT 120
  ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
endif()

set(_thread_args)
if(CHANNEL_DNS_MPI_ONLY)
  list(APPEND _thread_args --pure-mpi)
endif()
add_test(NAME thread_arguments
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_thread_args.py
    $<TARGET_FILE:channel_dns> ${CMAKE_SOURCE_DIR} ${_thread_args})
set_tests_properties(thread_arguments PROPERTIES TIMEOUT 120)
