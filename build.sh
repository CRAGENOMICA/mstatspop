#!/bin/bash
set -e

BUILD_TYPE=${BUILD_TYPE:-Release}
CFLAGS=${CFLAGS:-""}

# check if the build directory exists
mkdir -p build
cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_C_FLAGS="$CFLAGS" -DCMAKE_CXX_FLAGS="$CFLAGS" -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -S./ -B./build 
cmake --build ./build --config Release  -j 14 --