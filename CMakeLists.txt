cmake_minimum_required(VERSION 3.5)

project(BdfCGAL LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
    -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

find_package(CGAL COMPONENTS ImageIO)

set(SOURCES main.cpp nastran_export.hpp crap.hpp)

if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )

    find_package(Boost)

    # Activate concurrency ? (turned OFF by default)
    option(CGAL_ACTIVATE_CONCURRENT_MESH_3
        "Activate parallelism in Mesh_3"
        OFF)

    # And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
    if( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
        add_definitions( -DCGAL_CONCURRENT_MESH_3 )
        find_package( TBB REQUIRED )
    else( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
        option( LINK_WITH_TBB
            "Link with TBB anyway so we can use TBB timers for profiling"
            ON)
        if( LINK_WITH_TBB )
            find_package( TBB )
        endif( LINK_WITH_TBB )
    endif()

    # Use Eigen
    find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
    if (NOT EIGEN3_FOUND)
        message(STATUS "This project requires the Eigen library, and will not be compiled.")
        return()
    else()
        include( ${EIGEN3_USE_FILE} )
    endif()

    find_package(VTK QUIET COMPONENTS vtkImagingGeneral vtkIOImage NO_MODULE)
    if(VTK_FOUND)
        if(VTK_USE_FILE)
            include(${VTK_USE_FILE})
        endif()
        if ("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
            message(STATUS "VTK found")
            if(TARGET VTK::IOImage)
                set(VTK_LIBRARIES VTK::ImagingGeneral VTK::IOImage)
            endif()
        else()
            message(STATUS "VTK version 6.0 or greater is required")
        endif()
    else()
        message(STATUS "VTK was not found")
    endif()

    add_executable(BdfCGAL ${SOURCES})
    target_link_libraries(BdfCGAL CGAL::CGAL)

else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
