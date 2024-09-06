cmake_minimum_required(VERSION 3.12)

project(custom_matrix_types)

set(test_list "Tests/constexpr_basic_test.cpp" )

set(CMAKE_CXX_STANDARD 17)
add_executable(tests ${test_list})

#target_compile_features(tests PUBLIC cxx_std_17)