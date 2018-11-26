# ===============
# === Project ===
# ===============
cmake_minimum_required(VERSION 3.10)
project(ExaStamp VERSION 1.1 LANGUAGES CXX)

# ========================================
# === Compiler toolchain configuration ===
# ========================================
# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)


# ExaStamp binary 
add_executable(soatest SOAAccessor.cpp)

# Install commands
install(TARGETS soatest DESTINATION bin)

