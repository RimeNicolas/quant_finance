cmake_minimum_required(VERSION 3.0.0)
project(financial_engine VERSION 0.1.0)

include(CTest)
enable_testing()

file (GLOB financial_engine_SRC 
    "src/*.h"
    "src/*.cpp"
)

add_executable(financial_engine ${financial_engine_SRC})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
