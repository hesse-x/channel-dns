add_library(channel_dns_core STATIC
  core/params.cpp
  core/mesh.cpp
  core/mpi_context.cpp
)
target_include_directories(channel_dns_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_core PUBLIC MPI::MPI_CXX PkgConfig::FFTW)
if(CHANNEL_DNS_PURE_MPI)
  target_compile_definitions(channel_dns_core PUBLIC CHANNEL_DNS_PURE_MPI=1)
else()
  target_sources(channel_dns_core PRIVATE core/cpu_device.cpp)
  target_link_libraries(channel_dns_core PRIVATE Threads::Threads)
endif()

add_library(channel_dns_fft STATIC
  fft/transforms.cpp
  fft/wisdom.cpp
)
target_include_directories(channel_dns_fft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_fft PUBLIC channel_dns_core PkgConfig::FFTW)
if(NOT CHANNEL_DNS_PURE_MPI)
  target_link_libraries(channel_dns_fft PRIVATE ${FFTW_THREADS_LIBRARY})
endif()

add_library(channel_dns_solver STATIC
  solver/vector_ops.cpp
  solver/eqsolve_pv.cpp
)
target_include_directories(channel_dns_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_solver PUBLIC channel_dns_core channel_dns_fft)

add_library(channel_dns_timestep STATIC
  solver/timestepper.cpp
)
target_include_directories(channel_dns_timestep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_timestep PUBLIC channel_dns_solver channel_dns_fft)

add_library(channel_dns_io STATIC
  io/field_io.cpp
)
target_include_directories(channel_dns_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_io PUBLIC channel_dns_core channel_dns_fft channel_dns_solver)

add_library(channel_dns_sim STATIC
  sim/initial_bc.cpp
  sim/statistics.cpp
  sim/diagnostics.cpp
  sim/simulation.cpp
)
target_include_directories(channel_dns_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(channel_dns_sim PUBLIC
  channel_dns_core channel_dns_fft channel_dns_solver channel_dns_timestep channel_dns_io)

add_executable(channel_dns main.cpp)
target_link_libraries(channel_dns PRIVATE channel_dns_sim)
