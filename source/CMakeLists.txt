# cmake_minimum_required(VERSION 3.12)
cmake_minimum_required(VERSION 3.5)

# Project name
project(PSEOPT VERSION 0.1)

# C++14 is required to build this project
set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Position Independent Code (need for python)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()



###########################################
#Boost
############################################ 
option( PESOPT_WITH_BOOST "Boost" ON )
if ( PESOPT_WITH_BOOST )

    set(Boost_NO_BOOST_CMAKE ON )
    set(Boost_NO_SYSTEM_PATHS ON )

    #alternative: set BOOST_ROOT directory
     set(PESOPT_BOOST_DIR "" CACHE PATH "The directory where BOOST library is installed")
     if(PESOPT_BOOST_DIR)
         message( "set boost_Root to" ${PESOPT_BOOST_DIR} ) 
         set(BOOST_ROOT ${PESOPT_BOOST_DIR})
         #set(BOOST_INCLUDEDIR ${PESOPT_BOOST_DIR}/include/boost-1_75)
         set(BOOST_INCLUDEDIR ${PESOPT_BOOST_DIR}/include)
         set(BOOST_LIBRARYDIR ${PESOPT_BOOST_DIR}/lib)
     endif(PESOPT_BOOST_DIR)


    set (Boost_USE_STATIC_LIBS OFF) #must be OFF if one wants to build shared libraries
    #set (Boost_USE_STATIC_LIBS ON)
    set (Boost_USE_MULTITHREADED ON)  
    set (Boost_USE_STATIC_RUNTIME OFF) 
    
    set(BOOST_MIN_VERSION "1.75.0")
# 
#    find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
    find_package (Boost COMPONENTS system filesystem program_options log log_setup ) 
    #include_directories (${Boost_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    
    #    set_property(TARGET ${Boost_INCLUDE_DIRS} PROPERTY POSITION_INDEPENDENT_CODE ON)
    add_definitions(-DPESOPT_WITH_BOOST)
    
    ###
    message("-----------------------------------------")
    message("Use PESOPT with BOOST")
    message(STATUS "Boost version: ${Boost_VERSION}")
    message("-----------------------------------------")
endif ( PESOPT_WITH_BOOST )


###########################################
#Eigen
############################################ 
option( PESOPT_WITH_EIGEN "Eigen Linear Algebra Library" ON )
if( PESOPT_WITH_EIGEN )
    ### Alternative 1: directly use path to local version
    # include_directories ("${PROJECT_SOURCE_DIR}/externalLibs/eigen")
    ###
    ### Alternative 2:
    find_path(EIGEN_INCLUDE_DIR Eigen/Core PATHS "${PROJECT_SOURCE_DIR}/externalLibs/eigen" NO_DEFAULT_PATH)                              
    include_directories(${EIGEN_INCLUDE_DIR})
    ###
    # if(NOT EIGEN_INCLUDE_DIR)
    #   message(FATAL_ERROR "Could not find the Eigen library. Please make sure it is installed and included in the CMake prefix path.")
    # endif(NOT EIGEN_INCLUDE_DIR)
    add_definitions(-DPESOPT_WITH_EIGEN)
    
    ###
    message("-----------------------------------------")
    message("Use PESOPT with EIGEN")
    message("-----------------------------------------")
endif( PESOPT_WITH_EIGEN )
               
###########################################
#Ceres
###########################################
# option( PESOPT_WITH_CERES "Ceres Solver from Google" OFF )
# if( PESOPT_WITH_CERES )
#  find_package(Ceres REQUIRED)
#  include_directories(${CERES_INCLUDE_DIRS})
# #  add_definitions(-DPESOPT_WITH_CERES)
#  link_directories(${CERES_LIBRARY_DIR})
# endif( PESOPT_WITH_CERES)


###########################################
#IPOPT
###########################################
option( PESOPT_WITH_IPOPT "Interior Point Solver" ON )
if( PESOPT_WITH_IPOPT )

    find_path(IPOPT_INCLUDE_DIR IpTNLP.hpp PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt/build/include/coin" NO_DEFAULT_PATH )  
    find_library(IPOPT_LIBRARY ipopt PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt/build/lib" NO_DEFAULT_PATH ) 
    
    set(IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIR}" )
    set(IPOPT_LIBRARIES "${IPOPT_LIBRARY}")    
    add_definitions(-DPESOPT_WITH_IPOPT)
    #wrong config under Debian workaround    
    add_definitions( -DHAVE_CSTDDEF )    
    
    include(FindPackageHandleStandardArgs)   
    # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE if all listed variables are TRUE   
    find_package_handle_standard_args(IPOPT  DEFAULT_MSG  IPOPT_LIBRARY IPOPT_INCLUDE_DIR) 
    
    include_directories(${IPOPT_INCLUDE_DIR})
    link_directories(${IPOPT_LIBRARY})
    
    ###
    message("-----------------------------------------")
    message("Use PESOPT with IPOPT")
    message("-----------------------------------------")
endif( PESOPT_WITH_IPOPT )


###########################################
#OpenMP
############################################  
option( PESOPT_WITH_OPENMP "OpenMP for parallel computing" ON )
if ( PESOPT_WITH_OPENMP )
    find_package(OpenMP)
    add_definitions( -DPESOPT_WITH_OPENMP )
    ###
    message("-----------------------------------------")
    message("Use PESOPT with OPENMP")
    message("-----------------------------------------")
endif( PESOPT_WITH_OPENMP )



###########################################
#MPI
###########################################
option( PESOPT_WITH_MPI "MPI for parallel computing" OFF )
if ( PESOPT_WITH_MPI )
    find_package(MPI)
    find_path(MPI_INCLUDE_DIR mpi.h PATHS "${PROJECT_SOURCE_DIR}/externalLibs/openmpi-4.0.0/build/ompi/include/" NO_DEFAULT_PATH)                              
    include_directories(${MPI_INCLUDE_DIR})
    # #  include_directories(${MPI_INCLUDE_PATH})
    #  include_directories ("${PROJECT_SOURCE_DIR}/externalLibs/openmpi-4.0.0/build/ompi/include/")
    #  #for OPENMPI
    # #  add_definitions(-DOMPI_SKIP_MPICXX) 
    # #  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    # #  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    # #  include_directories(MPI_INCLUDE_PATH)
    add_definitions( -DPESOPT_WITH_MPI )
    ###
    message("-----------------------------------------")
    message("Use PESOPT with MPI")
    message("-----------------------------------------")
endif( PESOPT_WITH_MPI )

###########################################
#SuiteSparse
############################################ 
option( PESOPT_WITH_SUITESPARSE "SuiteSparse provides direct solvers for linear systems" ON )
if ( PESOPT_WITH_SUITESPARSE )
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
# find_package(SUITESPARSE REQUIRED)
# if(NOT SUITESPARSE_INCLUDE_DIRS)
#   message(FATAL_ERROR "Could not find the SUITESPARSE library. Please make sure it is installed and included in the CMake prefix path.")
# endif(NOT SUITESPARSE_INCLUDE_DIRS)

find_path(CHOLMOD_INCLUDE_DIR NAMES cholmod.h PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/include" NO_DEFAULT_PATH )
if ( CHOLMOD_INCLUDE_DIR )
  list ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
endif( CHOLMOD_INCLUDE_DIR )

find_path( SUITESPARSE_LIBRARY_DIR NAMES libcholmod.so libcholmod.a PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/lib" NO_DEFAULT_PATH )
if ( SUITESPARSE_LIBRARY_DIR )
  list ( APPEND SUITESPARSE_LIBRARIES amd)
  list ( APPEND SUITESPARSE_LIBRARIES btf)
  list ( APPEND SUITESPARSE_LIBRARIES camd)
  list ( APPEND SUITESPARSE_LIBRARIES ccolamd)
  list ( APPEND SUITESPARSE_LIBRARIES cholmod)
  list ( APPEND SUITESPARSE_LIBRARIES colamd)
  #list ( APPEND SUITESPARSE_LIBRARIES csparse)
  list ( APPEND SUITESPARSE_LIBRARIES cxsparse)
  list ( APPEND SUITESPARSE_LIBRARIES klu)
  #list ( APPEND SUITESPARSE_LIBRARIES spqr)
  list ( APPEND SUITESPARSE_LIBRARIES umfpack)

  # from version 4 on an additional config library is needed
#   find_library( SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig PATHS ${SUITESPARSE_LIBRARY_DIR})
  find_library( SUITESPARSE_CONFIG_LIB NAMES libsuitesparseconfig.so PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/lib" NO_DEFAULT_PATH )
  if ( EXISTS ${SUITESPARSE_CONFIG_LIB} )
    list ( APPEND SUITESPARSE_LIBRARIES ${SUITESPARSE_CONFIG_LIB} )
  endif ( EXISTS ${SUITESPARSE_CONFIG_LIB} )

  # Metis is optional
#   find_library( SUITESPARSE_METIS_LIBRARY NAMES metis PATHS ${SUITESPARSE_LIBRARY_DIR} )
#   if (SUITESPARSE_METIS_LIBRARY)                   
#     list ( APPEND SUITESPARSE_LIBRARIES metis)
#   endif(SUITESPARSE_METIS_LIBRARY)
  
#   find_library( SUITESPARSE_METIS_LIBRARY metis coinmetis PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt/build/lib/" NO_DEFAULT_PATH )
# #   find_library( SUITESPARSE_METIS_LIBRARY NAMES metis PATHS ${SUITESPARSE_LIBRARY_DIR} )
#   if (SUITESPARSE_METIS_LIBRARY)                   
#     list ( APPEND SUITESPARSE_LIBRARIES metis)
#   endif(SUITESPARSE_METIS_LIBRARY)
  
  #spqr is optional 
#   if(EXISTS  "${CHOLMOD_INCLUDE_DIR}/SuiteSparseQR.hpp")
#     set(SUITESPARSE_SPQR_VALID TRUE CACHE BOOL "SuiteSparseSPQR valid")
#   else()
#     set(SUITESPARSE_SPQR_VALID false CACHE BOOL "SuiteSparseSPQR valid")
#   endif()
else( SUITESPARSE_LIBRARY_DIR )  
  MESSAGE ( WARNING "SUITESPARSE_LIBRARY_DIR not found!" )
endif( SUITESPARSE_LIBRARY_DIR )  

if (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set(SUITESPARSE_FOUND TRUE)
else (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set( SUITESPARSE_FOUND FALSE )
endif (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

include_directories(${SUITESPARSE_INCLUDE_DIRS})
link_directories(${SUITESPARSE_LIBRARY_DIR})

add_definitions( -DPESOPT_WITH_SUITESPARSE )
endif( PESOPT_WITH_SUITESPARSE )

###########################################
#VTK
############################################  
option( PESOPT_WITH_VTK "VTK Library" ON )
if ( PESOPT_WITH_VTK )
    ### Alternative 1:
    #set(VTK_DIR "${PROJECT_SOURCE_DIR}/externalLibs/VTK/build" CACHE PATH "VTK directory override" FORCE)
    ###
    ### Alternative 2:
    find_path(VTK_DIR VTKConfig.cmake PATHS "${PROJECT_SOURCE_DIR}/externalLibs/VTK/build" NO_DEFAULT_PATH ) 
    ###
    find_package(VTK)
    if( ${VTK_MAJOR_VERSION} GREATER 8 )
      message( "VTK_MAJOR_VERSION bigger 8" )
    else ( )
      message( "VTK_MAJOR_VERSION smaller 9" )
      include(${VTK_USE_FILE})
    endif ( )
    add_definitions( -DPESOPT_WITH_VTK )
    ###
    message("-----------------------------------------")
    message("Use PESOPT with VTK:")
    message(VTK_MAJOR_VERSION = ${VTK_MAJOR_VERSION})
    message(VTK_MINOR_VERSION = ${VTK_MINOR_VERSION})
    message("-----------------------------------------")
endif( PESOPT_WITH_VTK )




############################################
#Include Projects from cmake.selection
############################################ 
if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )
  FILE ( COPY cmake.selection.default DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  FILE ( RENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection.default ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )
endif ( )

include ( ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )

foreach ( SELECTED_MODULE ${SELECTED_MODULES} )
 include_directories ( ${SELECTED_MODULE} )
endforeach ( SELECTED_MODULE )

foreach ( SELECTED_EXTERNAL ${SELECTED_EXTERNALS} )
 include_directories ( ${SELECTED_EXTERNAL} )
endforeach ( SELECTED_EXTERNAL )

foreach ( SELECTED_PROJECT ${SELECTED_PROJECTS} )
  add_subdirectory ( ${SELECTED_PROJECT} )
endforeach ( SELECTED_PROJECT )

foreach ( SELECTED_EXERCISE ${SELECTED_EXERCISES} )
  add_subdirectory ( ${SELECTED_EXERCISE} )
endforeach ( SELECTED_EXERCISE )

foreach ( SELECTED_TESTPROJECT ${SELECTED_TESTPROJECTS} )
  add_subdirectory ( ${SELECTED_TESTPROJECT} )
endforeach ( SELECTED_TESTPROJECT )

foreach ( SELECTED_TOOL ${SELECTED_TOOLS} )
  add_subdirectory ( ${SELECTED_TOOL} )
endforeach ( SELECTED_TOOL )



############################################
#automatically build documention with doxygen
############################################ 
# option ( USE_DOXYGEN "" OFF )
# if ( USE_DOXYGEN )
#   find_package ( Doxygen )
#   find_package ( LATEX )
#   if ( DOXYGEN_FOUND AND PDFLATEX_COMPILER AND PYTHONINTERP_FOUND )
#     #! \cmakeoption{Make Doxygen use MathJax to display Latex formulas instead of generating PNGs. Looks much better, in particular, on high-DPI displays.,OFF}
#     option ( USE_MATHJAX_FOR_DOXYGEN "" OFF )
#       #! \cmakeoption{Relative (to doc/lib, doc/all) path of MathJax,""}
#       set ( DOXYGEN_MATHJAX_RELPATH "" CACHE STRING "relative (to doc/lib, doc/all) path of mathjax")
#     if ( USE_MATHJAX_FOR_DOXYGEN )
#       if ( DOXYGEN_MATHJAX_RELPATH STREQUAL "" )
# 	message ( FATAL_ERROR "USE_MATHJAX_FOR_DOXYGEN is ON, but no DOXYGEN_MATHJAX_RELPATH has been set!" )
#       endif ( DOXYGEN_MATHJAX_RELPATH STREQUAL "" )
#       set( DOXYGEN_USE_MATHJAX "YES" )
#     else ( )
#       set( DOXYGEN_USE_MATHJAX "NO" )
#     endif ( )
#     # doclib
#     FILE ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ )
#     if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc )
#       FILE ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ )
#     endif ( )
#     make_directory ( ${CMAKE_CURRENT_BINARY_DIR}/doc/lib )
#     # input files to scan for documentation
#     foreach ( SELECTED_MODULE ${SELECTED_MODULES} )
#       # add module directory to the input file list
#       set ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_MODULE} )
#     endforeach ( SELECTED_MODULE )
#     # add doc folder to the input file list
#     set ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc )
#     # if the cmake doxygen filter is available, use it for CMakeLists.txt files
#     if ( USE_CMAKE_DOXYGEN_FILTER )
#       set ( DOXYGEN_FILTERS "*.txt=${CMAKE_CURRENT_BINARY_DIR}/internal/external/cmakeDoxygenFilter/CMakeDoxygenFilter" )
#       set ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt )
#     endif ( USE_CMAKE_DOXYGEN_FILTER )
#     STRING( REGEX REPLACE ";" " " DOXYGEN_INPUT_LIST "${DOXYGEN_INPUT}" )
#     # preprocessor defines
#     get_property( DEFS DIRECTORY PROPERTY COMPILE_DEFINITIONS )
#     STRING( REGEX REPLACE ";" " " DOXYGEN_PREDEFINED "${DEFS}" )
#     # other stuff for doxygen
#     set( DOXYGEN_ENABLED_SECTIONS "DOC_LIB_ONLY" )
#     set( DOXYGEN_HTML_OUTPUT "lib" )
#     set( DOXYGEN_MAN_OUTPUT "lib/man" )
#     # prepare doxygen configuration file
#     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.cmake ${CMAKE_CURRENT_BINARY_DIR}/doc/lib/doxygen.conf )
#     # add target
#     ADD_CUSTOM_TARGET ( doclib
#                                 ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${DOXYGEN_EXECUTABLE} lib/doxygen.conf
#                         WORKING_DIRECTORY doc
#     )
# 
#     # docall
#     make_directory ( ${CMAKE_CURRENT_BINARY_DIR}/doc/all )
#     # additional input files to scan for documentation
#     FOREACH ( SELECTED_PROJECT ${SELECTED_PROJECTS} ${SELECTED_STDPROJECTS} )
#       # add project directory to the input file list
#       set ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_PROJECT} )
#     ENDFOREACH ( SELECTED_PROJECT )
#     STRING( REGEX REPLACE ";" " " DOXYGEN_INPUT_LIST "${DOXYGEN_INPUT}" )
#     # other stuff for doxygen
#     set( DOXYGEN_ENABLED_SECTIONS "" )
#     set( DOXYGEN_HTML_OUTPUT "all" )
#     set( DOXYGEN_MAN_OUTPUT "all/man" )
#     # prepare doxygen configuration file
#     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.cmake ${CMAKE_CURRENT_BINARY_DIR}/doc/all/doxygen.conf )
#     # add target
#     ADD_CUSTOM_TARGET ( docall
#                                 ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  manual.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
#                         COMMAND ${DOXYGEN_EXECUTABLE} all/doxygen.conf
#                         WORKING_DIRECTORY doc
#     )
# 
#     # if the cmake doxygen filter is enabled it has to build before the doc* targets are executed
#     if ( USE_CMAKE_DOXYGEN_FILTER )
#       add_dependencies( doclib CMakeDoxygenFilter )
#       add_dependencies( docall CMakeDoxygenFilter )
#     endif ()
# 
#     # Table of projects generation is only available for linux
#     ADD_CUSTOM_COMMAND( TARGET docall POST_BUILD COMMAND bash ARGS ${CMAKE_CURRENT_SOURCE_DIR}/util/makeTableOfProjects.sh ${SELECTED_PROJECTS} ${SELECTED_STDPROJECTS} )
# 
#     # cleanup folder doc on make clean
#     set_PROPERTY( DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES doc )
#   else ( )
#     set ( WARNING_TEXT "" )
#     if ( NOT DOXYGEN_FOUND  )
#       set ( WARNING_TEXT "${WARNING_TEXT}Doxygen not available, " )
#     endif ( )
#     if ( NOT PDFLATEX_COMPILER )
#       set ( WARNING_TEXT "${WARNING_TEXT}pdflatex compiler not available, " )
#     endif ( )
#     if ( NOT PYTHONINTERP_FOUND )
#       set ( WARNING_TEXT "${WARNING_TEXT}python interpreter not available, " )
#     endif ( )
#     MESSAGE ( WARNING "${WARNING_TEXT}not generating doc targets." )
#   endif ( )
# endif ( )




#  To let the script know where the source and thus the error parse scripts are sneak in the path as argument.
set_property ( GLOBAL PROPERTY RULE_LAUNCH_COMPILE "bash ${CMAKE_CURRENT_SOURCE_DIR}/util/cmakeParseError.sh ${CMAKE_CURRENT_SOURCE_DIR}" )
