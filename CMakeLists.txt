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

