cmake_minimum_required(VERSION 4.0)
project(Toolbox)

set(CMAKE_CXX_STANDARD 20)

set(TOOLBOX_SOURCES
        baseline/BinBaseline.cpp
        baseline/BinBaseline.h
)

add_executable(Toolbox
        ${TOOLBOX_SOURCES}
        main.cpp
        var_needed.h
        models/Bin.h
)

add_library(ToolboxDLL
        ${TOOLBOX_SOURCES}
        api/API.cpp
        api/API.h
)
