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

find_package(OpenMP REQUIRED)

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

macro(GenerateBenchmark A C)
  add_executable(soatlbenchmark_${A}_${C} tests/benchmark.cpp)
  target_include_directories(soatlbenchmark_${A}_${C} PUBLIC include)
  target_compile_options(soatlbenchmark_${A}_${C} PUBLIC ${OpenMP_CXX_FLAGS} -DTEST_ALIGNMENT=${A} -DTEST_CHUNK_SIZE=${C} )
  target_link_libraries(soatlbenchmark_${A}_${C} ${OpenMP_LIBRARIES})
endmacro()

GenerateBenchmark(64 16)
GenerateBenchmark(64 8)
GenerateBenchmark(64 4)
GenerateBenchmark(64 2)
GenerateBenchmark(64 1)
GenerateBenchmark(32 8)
GenerateBenchmark(32 4)
GenerateBenchmark(32 2)
GenerateBenchmark(32 1)
GenerateBenchmark(16 4)
GenerateBenchmark(16 2)
GenerateBenchmark(16 1)
GenerateBenchmark(8 1)
GenerateBenchmark(1 1)

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

