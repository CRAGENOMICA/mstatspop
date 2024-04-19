cmake_minimum_required(VERSION 3.5.0)
project(mstatspop VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()


add_compile_options("-DLOG_USE_COLOR")
IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
add_compile_options("-DDEBUG")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

file (GLOB_RECURSE SOURCE_FILE CONFIGURE_DEPENDS  "mstatspop/*.c" "lib/*.c")

add_executable(mstatspop ${SOURCE_FILE})
target_include_directories(mstatspop PRIVATE "lib")
target_compile_options(mstatspop PRIVATE -DinGSL=1)
target_link_libraries(mstatspop  -lm -lgsl -lgslcblas -lz)

# Assuming HTSlib was installed in /usr/local
find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS /usr/local/include)
find_library(HTSLIB_LIBRARY NAMES hts PATHS /usr/local/lib)

## index_create 
add_executable(create_index index/main.c)
target_include_directories(create_index PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(create_index PRIVATE ${HTSLIB_LIBRARY})



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
