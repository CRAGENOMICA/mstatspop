cmake_minimum_required(VERSION 3.5.0)
project(mstatspop VERSION 1.1.0 LANGUAGES C CXX)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

# include(CTest)
# enable_testing()



# Get the current date and time
string(TIMESTAMP CURRENT_DATE "%Y%m%d")
string(TIMESTAMP CURRENT_TIME "%H%M%S")


# Concatenate date and time to create a build number
set(BUILD_NUMBER "${CURRENT_DATE}${CURRENT_TIME}")

# Define the version number
set(VERSION_NUMBER "${mstatspop_VERSION}")

# Define preprocessor macros for date, time, and build number
add_definitions(-DCURRENT_DATE=\"${CURRENT_DATE}\")
add_definitions(-DCURRENT_TIME=\"${CURRENT_TIME}\")
add_definitions(-DBUILD_NUMBER=\"${BUILD_NUMBER}\")
add_definitions(-DVERSION_NUMBER=\"${VERSION_NUMBER}\")

# Write version number to a text file
file(WRITE "${CMAKE_BINARY_DIR}/version.txt" "${VERSION_NUMBER}")





add_compile_options("-DLOG_USE_COLOR")
IF(CMAKE_BUILD_TYPE MATCHES Debug)

# Set compiler options for warnings and AddressSanitizer
# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     # add_compile_options(-Wall -Wextra -pedantic -fsanitize=address -g)
#     add_compile_options(-Wall -Wextra -pedantic -fsanitize=address -g -O0)
#     add_link_options(-fsanitize=address)
# endif()

add_compile_options("-DDEBUG")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

file (GLOB_RECURSE SOURCE_FILE CONFIGURE_DEPENDS  "mstatspop/*.c")
## explude mstatspop/tfa_index.c and mstatspop/tfa_convert.c
list(REMOVE_ITEM SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mstatspop/tfa_index.c")
list(REMOVE_ITEM SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mstatspop/tfa_merge.c")
list(REMOVE_ITEM SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mstatspop/tfa_convert.c")

# log library
add_library(logging "lib/log.c")


#add_library(zutil "lib/zutil.c")


# Assuming HTSlib was installed in /usr/local
find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS /usr/local/include)
find_library(HTSLIB_LIBRARY NAMES hts PATHS /usr/local/lib)

## find htslib on windows NEED TO BE TESTED
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "C:/Program Files/htslib/include")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "C:/Program Files/htslib/lib")
endif()
## find htslib on macos NEED TO BE TESTED
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/usr/local/Cellar/htslib")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/usr/local/Cellar/htslib")
endif()
# /opt/homebrew/opt/htslib
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/opt/homebrew/opt/htslib/include")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/opt/homebrew/opt/htslib/lib")
endif()
## if not search in  /opt/homebrew/Cellar/htslib/1.20
if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib/sam.h HINTS "/opt/homebrew/Cellar/htslib")
  find_library(HTSLIB_LIBRARY NAMES hts PATHS "/opt/homebrew/Cellar/htslib") 
endif()

if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
  message(FATAL_ERROR "Could not find htslib")
endif()

## set GLS_INCLUDE_DIR
# find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h HINTS /usr/local/include)
# find_library(GSL_LIBRARY NAMES gsl PATHS /usr/local/lib)

# if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY)
#   find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h HINTS "C:/Program Files/gsl/include")
#   find_library(GSL_LIBRARY NAMES gsl PATHS "C:/Program Files/gsl/lib")
# endif()
# ## find gls on macos NEED TO BE TESTED
# if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY)
#   find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h HINTS "/usr/local/Cellar/gsl")
#   find_library(GSL_LIBRARY NAMES gsl PATHS "/usr/local/Cellar/gsl")
# endif()
# # /opt/homebrew/opt/gsl
# if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY)
#   find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h HINTS "/opt/homebrew/opt/gsl/include")
#   find_library(GSL_LIBRARY NAMES gsl PATHS "/opt/homebrew/opt/gsl/lib")
# endif()
# ## if not search in  /opt/homebrew/Cellar/gsl
# if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY)
#   find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h HINTS "/opt/homebrew/Cellar/gsl")
#   find_library(GSL_LIBRARY NAMES gsl PATHS "/opt/homebrew/Cellar/gsl")
# endif()

# if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY)
#   message(FATAL_ERROR "Could not find gsl")
# endif()

FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

include_directories(${HTSLIB_INCLUDE_DIR} ${GLS_INCLUDE_DIR} ${GSLCBLAS_INCLUDE_DIRS})

add_executable(mstatspop ${SOURCE_FILE})
target_include_directories(mstatspop PRIVATE "lib")
target_compile_options(mstatspop PRIVATE -DinGSL=1)
target_link_libraries(mstatspop  ${HTSLIB_LIBRARY}  ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES} -lm  -lz logging )
# target_link_libraries(mstatspop PRIVATE  )



## tfa_index 
add_executable(tfa_index mstatspop/tfa_index.c mstatspop/tfasta.c mstatspop/files_util.c)
target_include_directories(tfa_index PRIVATE "lib")
target_include_directories(tfa_index PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(tfa_index PRIVATE ${HTSLIB_LIBRARY} logging)

## tfa_merge
add_executable(tfa_merge mstatspop/tfa_merge.c mstatspop/tfasta.c mstatspop/files_util.c)
target_include_directories(tfa_merge PRIVATE "lib")
target_include_directories(tfa_merge PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(tfa_merge PRIVATE ${HTSLIB_LIBRARY} logging)

## ms
add_executable(ms sources_msHudson/ms.c  sources_msHudson/streec.c sources_msHudson/rand2.c)
target_include_directories(ms PRIVATE sources_msHudson )
target_link_libraries(ms -lm )

## gVCF2tFasta
#file (GLOB_RECURSE SOURCE_FILE CONFIGURE_DEPENDS  "gvcf2tfasta/*.cpp")

#add_executable(gVCF2tFasta ${SOURCE_FILE} )
#target_include_directories(gVCF2tFasta PRIVATE "lib")
#target_link_libraries(gVCF2tFasta PRIVATE logging ${HTSLIB_LIBRARY})


## tfa_fadix 
#add_executable(tfa_fadix ctests/tfa_fadix.c)
#target_include_directories(tfa_fadix PRIVATE "lib")
#target_include_directories(tfa_fadix PRIVATE ${HTSLIB_INCLUDE_DIR})
#target_link_libraries(tfa_fadix PRIVATE ${HTSLIB_LIBRARY} logging)

## tfa_test
#add_executable(test_tfasta ctests/test_tfasta.c mstatspop/zutil.c mstatspop/zindex.c mstatspop/get_tfadata.c mstatspop/usegff.c mstatspop/get_obsdatastats.c)
#target_include_directories(test_tfasta PRIVATE "lib")
#target_include_directories(test_tfasta PRIVATE "mstatspop")
#target_include_directories(test_tfasta PRIVATE ${HTSLIB_INCLUDE_DIR})
#target_link_libraries(test_tfasta PRIVATE -lm ${HTSLIB_LIBRARY} -lz logging)

## tfa_convert
#add_executable(tfa_convert ctests/tfa_convert.c)
#target_include_directories(tfa_convert PRIVATE "lib")
#target_include_directories(tfa_convert PRIVATE ${HTSLIB_INCLUDE_DIR})
#target_link_libraries(tfa_convert PRIVATE ${HTSLIB_LIBRARY} logging)

## test_zutil
#add_executable(test_zutil ctests/test_zutil.c mstatspop/zutil.c mstatspop/zindex.c)
## add the include directories to the test_zutil
#target_include_directories(test_zutil PRIVATE "mstatspop")
#target_link_libraries(test_zutil -lz )

## test_compress
#add_executable(test_compress ctests/test_compress.c mstatspop/zutil.c mstatspop/zindex.c)
## add the include directories to the test_compress
#target_include_directories(test_compress PRIVATE "mstatspop")
#target_link_libraries(test_compress -lz )


## test_compress
#add_executable(gun ctests/gun.c)
## add the include directories to the test_compress
#target_link_libraries(gun -lz )

## test_compress
#add_executable(gun1 ctests/gun1.c)
## add the include directories to the test_compress
# target_link_libraries(gun1 -lz )

## test_compress
#add_executable(zpipe ctests/zpipe.c)
## add the include directories to the test_compress
#target_link_libraries(zpipe -lz )


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
install(TARGETS mstatspop DESTINATION bin)
install(TARGETS ms DESTINATION bin)
install(TARGETS tfa_index DESTINATION bin)
install(TARGETS tfa_merge DESTINATION bin)
install(FILES README.md DESTINATION .)
include(CPack)
