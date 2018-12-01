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

# test binary 
add_executable(soatest tests/soatest.cpp)
target_include_directories(soatest PUBLIC include)
target_compile_options(soatest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatest ${OpenMP_LIBRARIES})

enable_testing()
add_test(NAME soatl_test1 COMMAND soatest 1000 0)
add_test(NAME soatl_test2 COMMAND soatest 1000 34523452)
add_test(NAME soatl_test3 COMMAND soatest 1000 1976)
add_test(NAME soatl_test4 COMMAND soatest 1000 234234234)

# Install commands
install(TARGETS soatest DESTINATION bin)

