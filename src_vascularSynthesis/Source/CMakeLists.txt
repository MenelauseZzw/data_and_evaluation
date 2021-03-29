cmake_minimum_required(VERSION 3.5)
project(VascuSynth)

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
endif()

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(VascuSynth VascuSynth.cpp SupplyMap.cpp OxygenationMap.cpp NodeTable.cpp VascularTree.cpp TreeDrawer.cpp)
target_link_libraries(VascuSynth ${ITK_LIBRARIES})
