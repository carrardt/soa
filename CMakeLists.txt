# ===============
# === Project ===
# ===============
cmake_minimum_required(VERSION 3.10)
project(ExaStamp VERSION 1.1 LANGUAGES CXX)

# ========================================
# === Compiler toolchain configuration ===
# ========================================
# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)

# packages and tools
find_package(OpenMP REQUIRED)
get_filename_component(BINUTILS_DIR ${CMAKE_LINKER} DIRECTORY)
find_file(SOATL_OBJDUMP objdump HINTS ${BINUTILS_DIR})

# Install commands
install(DIRECTORY include/soatl DESTINATION include)

# build tests binaries
add_executable(soatltest tests/soatest.cpp)
target_include_directories(soatltest PUBLIC include)
target_compile_options(soatltest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatltest ${OpenMP_CXX_LIB_NAMES})

add_executable(soatlcomputetest tests/computetest.cpp)
target_include_directories(soatlcomputetest PUBLIC include)
target_compile_options(soatlcomputetest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatlcomputetest ${OpenMP_CXX_LIB_NAMES})

add_executable(soatlserializetest tests/serializetest.cpp)
target_include_directories(soatlserializetest PUBLIC include)
target_compile_options(soatlserializetest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatlserializetest ${OpenMP_CXX_LIB_NAMES})

# register tests
enable_testing()
add_test(NAME soatl_test1 COMMAND soatltest 1000 0)
add_test(NAME soatl_test2 COMMAND soatltest 1000 34523452)
add_test(NAME soatl_test3 COMMAND soatltest 1000 1976)
add_test(NAME soatl_test4 COMMAND soatltest 1000 234234234)
add_test(NAME soatl_compute1 COMMAND soatlcomputetest 1000 0)
add_test(NAME soatl_compute2 COMMAND soatlcomputetest 1000 34523452)
add_test(NAME soatl_compute3 COMMAND soatlcomputetest 1000 1976)
add_test(NAME soatl_compute4 COMMAND soatlcomputetest 1000 234234234)
add_test(NAME soatl_serialize COMMAND soatlserializetest 10000)

# benchmarking
if(SOATL_OBJDUMP)
  add_custom_target(vecreport)
endif()

macro(GenerateBenchmark A C DPS SIMD OMPTOGGLE)
  if("${DPS}" STREQUAL "d")
    set(DP 1)
  else()
    set(DP 0)
  endif()
  if(${SIMD})
    set(VEC 1)
    set(VECS "_vec")
  else()
    set(VEC 0)
    set(VECS "")
  endif()
  if(${OMPTOGGLE})
    set(OMPVAL 1)
    set(OMPS "_omp")
  else()
    set(OMPS "")
    set(OMPVAL 0)
  endif()
  set(SUFFIX ${A}_${C}_${DPS}${VECS}${OMPS})
  add_executable(soatlbenchmark_${SUFFIX} tests/benchmark.cpp)
  target_include_directories(soatlbenchmark_${SUFFIX} PUBLIC include)
  target_compile_options(soatlbenchmark_${SUFFIX} PUBLIC ${OpenMP_CXX_FLAGS})
  target_compile_definitions(soatlbenchmark_${SUFFIX} PUBLIC -DTEST_USE_SIMD=${VEC} -DTEST_ALIGNMENT=${A} -DTEST_CHUNK_SIZE=${C} -DTEST_DOUBLE_PRECISION=${DP} -DTEST_ENABLE_OPENMP=${OMPVAL})
  target_link_libraries(soatlbenchmark_${SUFFIX} ${OpenMP_CXX_LIB_NAMES})
  # add perf tests
  add_test(NAME soatlbenchmark_${SUFFIX}_fa COMMAND soatlbenchmark_${SUFFIX} fa 10000000)
  add_test(NAME soatlbenchmark_${SUFFIX}_pfa COMMAND soatlbenchmark_${SUFFIX} pfa 10000000)
  # assembly analysis
  if(SOATL_OBJDUMP)
    add_custom_target(vecreport_${SUFFIX}
                  COMMAND ${CMAKE_COMMAND} -DBINARY_FILE="$<TARGET_FILE:soatlbenchmark_${SUFFIX}>" -DSOATL_OBJDUMP="${SOATL_OBJDUMP}"
                  -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/vecreport.cmake
                  DEPENDS soatlbenchmark_${SUFFIX})
    add_dependencies(vecreport vecreport_${SUFFIX})
  endif()
endmacro()

# benchmarks
GenerateBenchmark(64 16 d ON OFF)
GenerateBenchmark(64 16 d ON ON)
GenerateBenchmark(64 16 f ON OFF)
GenerateBenchmark(64 16 f ON ON)
GenerateBenchmark(32 8 d ON OFF)
GenerateBenchmark(32 8 f ON OFF)
GenerateBenchmark(16 4 d ON OFF)
GenerateBenchmark(16 4 f ON OFF)
GenerateBenchmark(1 1 d OFF OFF)
GenerateBenchmark(1 1 d OFF ON)
GenerateBenchmark(1 1 f OFF OFF)
GenerateBenchmark(1 1 f OFF ON)

