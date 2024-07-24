#!/bin/bash
set -e
# check if the build directory exists
mkdir -p build
cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -S./ -B./build 
cmake --build ./build --config Release  -j 14 --