# ============================================================================
# CMakeListst.txt file for MSTGRAPH PROJECT library                                        
# ============================================================================
cmake_minimum_required(VERSION 2.6)

# project name
project (TRACEID)
# find our custom CMake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TRACEID_SOURCE_DIR}/Modules)
#set (CMAKE_BUILD_TYPE Debug)
###################################
#SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_MULTITHREAD ON)
#SET( Boost_ADDITIONAL_VERSIONS "1.41.0" "1.40" )
find_package ( Boost  "1.41" COMPONENTS program_options thread )
################ TEST BOOST DIR 
MESSAGE (STATUS "BOOST FOUND IN " ${Boost_INCLUDE_DIRS})
MESSAGE (STATUS "BOOST LIB FOUND IN " ${Boost_LIBRARY_DIRS})
MESSAGE (STATUS "BOOST I WILL INCLUDE LIB  " ${Boost_LIBRARIES})
# contains the full path to the top level directory of your build tree 

##############################################
#find_package (GSL)
#find_package (GTS)
#find_package(ANN REQUIRED)
#MESSAGE (STATUS "ANN FOUND IN " ${ANN_INCLUDE_DIRS})
#MESSAGE (STATUS "ANN LIB FOUND IN " ${ANN_LIBRARIES})
#MESSAGE (STATUS "GTS FOUND IN " ${GTS_INCLUDE_DIR})
#MESSAGE (STATUS "Here GTS has: " ${CMAKE_GTS_CXX_FLAGS})
##############################################

MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )

# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command 
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
# Find all the sources for the LIB
#FILE(GLOB LIBSOURCES lib/*.cpp)

# Destination path for the lib
#SET(LIBRARY_OUTPUT_PATH lib)


if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "-Wno-deprecated -O3 -fopenmp ")
else(CMAKE_COMPILER_IS_GNUCXX)
        SET(CMAKE_CXX_FLAGS "-O3 -xW -Zp16 -ip -ansi_alias -ansi  -mp -vec-report3 -ipo " )	
endif(CMAKE_COMPILER_IS_GNUCXX)


#################
# Fixing missing glib2
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GTS_CXX_FLAGS}")

#################

MESSAGE(STATUS "BUILD TYPE: " ${CMAKE_BUILD_TYPE} )

# create library "armparallel"
#add_library (arman-tools SHARED ${LIBSOURCES})
#
# ------------------
# PROGRAMS          
# ------------------

# Executbale output path
set(EXECUTABLE_OUTPUT_PATH bin)
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

# Make sure the compiler can find include files from our Hello library.
include_directories (${TRACEID_SOURCE_DIR} ${TRACEID_SOURCE_DIR} ${Boost_INCLUDE_DIRS} )

# Make sure the linker can find the Hello library once it is built.
link_directories (${TRACEID_BINARY_DIR}  ${Boost_LIBRARY_DIRS} )


# ----------------------------------------------------------
# Make sure the compiler can find include files from our Hello library.

  add_executable (trace_id.x data_readers.cpp  options.cpp  trace_ids.cpp )

  # Link the executable to the Hello library.
  target_link_libraries (trace_id.x ${Boost_LIBRARIES} )


# ----------------------------------------------------------
# install target                                            

set (CMAKE_INSTALL_PREFIX $ENV{LOCAL})

# lib
#INSTALL(FILES  DESTINATION lib)

