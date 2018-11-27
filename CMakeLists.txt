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

# ExaStamp binary 
add_executable(soatest soatest.cpp)
target_compile_options(soatest PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(soatest ${OpenMP_LIBRARIES})


# Install commands
install(TARGETS soatest DESTINATION bin)

