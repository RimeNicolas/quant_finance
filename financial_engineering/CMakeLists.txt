cmake_minimum_required(VERSION 3.0.0)
project(financial_engine VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(financial_engine src/main.cpp src/hello.h src/hello.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
