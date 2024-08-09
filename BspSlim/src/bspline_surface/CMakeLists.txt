cmake_minimum_required(VERSION 3.0.2)
project(bspline_surface)


find_package(catkin REQUIRED COMPONENTS
  pcl_conversions
  pcl_msgs
  pcl_ros
  rospy
  std_msgs
)

find_package(PCL 1.10 REQUIRED)

# for Debug
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE "Release")


catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ptc_demo
 CATKIN_DEPENDS pcl_msgs pcl_ros roscpp rospy std_msgs pcl_conversions
#  DEPENDS system_lib
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)

## Declare a C++ library
# add_library(bsplinefitkd SHARED
#   src/bspline.cpp src/bspline_api.cpp
# )

link_directories(${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/lib)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(bsp_test src/bsp_test.cpp src/bspline_api.cpp src/bspline.cpp)
# add_executable(bsp_test src/bsp_test.cpp)


target_link_libraries(bsp_test
  # bsplinefit
  # bsplinefitkd
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}  
)

#############
## Install ##
#############