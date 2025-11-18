cmake_minimum_required(VERSION 4.0)
project(Toolbox)

set(CMAKE_CXX_STANDARD 20)

set(TOOLBOX_SOURCES
        src/baseline/BinBaseline.cpp
        src/baseline/BinBaseline.h
        src/baseline/AlsBaseline.cpp
        src/baseline/AlsBaseline.h
        src/models/elements/Bin.h
        src/models/ctx/BibSearchCtx.h
        src/models/cal/EC.h
        src/models/elements/FitOut.h
        src/models/cal/FWHMC.h
        src/models/ctx/LinearBGCtx.h
        src/models/elements/Nuclid.h
        src/nuclides/NuclidLibrary.h
        src/nuclides/NuclidLibrary.cpp
        src/finding/GammaMBibSearch.h
        src/finding/GammaMBibSearch.cpp
        #src/driver/GammaMBibSearchDriver.cpp
        src/driver/GammaMBibSearchTest.cpp
        src/finding/SecondDerivateSearch.cpp
        src/finding/SecondDerivateSearch.h
)

add_executable(Toolbox
        ${TOOLBOX_SOURCES}
        src/main.cpp
        src/driver/GammaMBibSearchTest.h
        src/driver/SecondDSearchTest.cpp
        src/driver/SecondDSearchTest.h
)

add_library(ToolboxDLL
        ${TOOLBOX_SOURCES}
        src/api/API.cpp
        src/api/API.h
)
