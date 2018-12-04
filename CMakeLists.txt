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

# base test : alignment, chunk boundaries, aliasing, resize, etc.
add_executable(soatltest tests/soatest.cpp)
target_include_directories(soatltest PUBLIC include)
target_compile_options(soatltest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatltest ${OpenMP_LIBRARIES})

add_executable(soatlcomputetest tests/computetest.cpp)
target_include_directories(soatlcomputetest PUBLIC include)
target_compile_options(soatlcomputetest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatlcomputetest ${OpenMP_LIBRARIES})

add_executable(soatlserializetest tests/serializetest.cpp)
target_include_directories(soatlserializetest PUBLIC include)
target_compile_options(soatlserializetest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatlserializetest ${OpenMP_LIBRARIES})

if(SOATL_OBJDUMP)
  add_custom_target(vecreport)
endif()

macro(GenerateBenchmark A C DPS SIMD)
  if("${DPS}" STREQUAL "d")
    set(DP 1)
  else()
    set(DP 0)
  endif()
  if(${SIMD})
    set(COMPFLAGS -DTEST_USE_SIMD=1)
    set(SUFFIX "vec")
  else()
    set(COMPFLAGS -DTEST_USE_SIMD=0)
    set(SUFFIX "scal")
  endif()
  add_executable(soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX} tests/benchmark.cpp)
  target_include_directories(soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX} PUBLIC include)
  target_compile_options(soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX} PUBLIC ${OpenMP_CXX_FLAGS} ${COMPFLAGS} -DTEST_ALIGNMENT=${A} -DTEST_CHUNK_SIZE=${C} -DTEST_DOUBLE_PRECISION=${DP})
  target_link_libraries(soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX} ${OpenMP_LIBRARIES})  
  if(SOATL_OBJDUMP)
    add_custom_target(vecreport_${A}_${C}_${DPS}_${SUFFIX}
                  COMMAND ${CMAKE_COMMAND} -DBINARY_FILE="$<TARGET_FILE:soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX}>" -DSOATL_OBJDUMP="${SOATL_OBJDUMP}" -P ${CMAKE_CURRENT_SOURCE_DIR}/vecreport.cmake
                  DEPENDS soatlbenchmark_${A}_${C}_${DPS}_${SUFFIX})
    add_dependencies(vecreport vecreport_${A}_${C}_${DPS}_${SUFFIX})
  endif()
endmacro()

# benchmarks
GenerateBenchmark(64 16 d ON)
GenerateBenchmark(64 16 f ON)
GenerateBenchmark(32 8 d ON)
GenerateBenchmark(32 8 f ON)
GenerateBenchmark(16 4 d ON)
GenerateBenchmark(16 4 f ON)
GenerateBenchmark(1 1 d OFF)
GenerateBenchmark(1 1 f OFF)

# tests
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
# Install commands
install(DIRECTORY include/soatl DESTINATION include)

