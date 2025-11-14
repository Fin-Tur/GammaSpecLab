cmake_minimum_required(VERSION 4.0)
project(Toolbox)

set(CMAKE_CXX_STANDARD 20)

set(TOOLBOX_SOURCES
        baseline/BinBaseline.cpp
        baseline/BinBaseline.h
        baseline/AlsBaseline.cpp
        baseline/AlsBaseline.h
        models/Bin.h
        models/BibSearchCtx.h
        models/EC.h
        models/FitOut.h
        models/FWHMC.h
        models/LinearBGCtx.h
        models/Nuclid.h
        nuclides/NuclidLibrary.h
        nuclides/NuclidLibrary.cpp
        finding/GammaMBibSearch.h
        finding/GammaMBibSearch.cpp
        driver/GammaMBibSearchDriver.cpp
        driver/GammaMBibSearchTest.cpp
        
)

add_executable(Toolbox
        ${TOOLBOX_SOURCES}
        main.cpp
        models/Bin.h
)

add_library(ToolboxDLL
        ${TOOLBOX_SOURCES}
        api/API.cpp
        api/API.h
)
