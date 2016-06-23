cmake_minimum_required(VERSION 2.8.3)
project(graph_slam_uk)
set(PCL_DIR "/usr/local/share/pcl-1.7/PCLConfig.cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules/")
## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
    roscpp
    roslaunch
    sensor_msgs
    tf
    laser_geometry
    message_filters
    nav_msgs
    pcl_ros
    geometry_msgs
    dynamic_slam_utils
    ndt_scanmatching2d
    tf_conversions
)

find_package(PCL 1.7 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options system)
find_package(Eigen3 REQUIRED)
# g2o #######################
find_package(G2O REQUIRED)
find_package(OpenMP)
## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)


## Uncomment this if the package has a setup.py. This macro ensures
## modules and global scripts declared therein get installed
## See http://ros.org/doc/api/catkin/html/user_guide/setup_dot_py.html
# catkin_python_setup()

################################################
## Declare ROS messages, services and actions ##
################################################

## To declare and build messages, services or actions from within this
## package, follow these steps:
## * Let MSG_DEP_SET be the set of packages whose message types you use in
##   your messages/services/actions (e.g. std_msgs, actionlib_msgs, ...).
## * In the file package.xml:
##   * add a build_depend and a run_depend tag for each package in MSG_DEP_SET
##   * If MSG_DEP_SET isn't empty the following dependencies might have been
##     pulled in transitively but can be declared for certainty nonetheless:
##     * add a build_depend tag for "message_generation"
##     * add a run_depend tag for "message_runtime"
## * In this file (CMakeLists.txt):
##   * add "message_generation" and every package in MSG_DEP_SET to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * add "message_runtime" and every package in MSG_DEP_SET to
##     catkin_package(CATKIN_DEPENDS ...)
##   * uncomment the add_*_files sections below as needed
##     and list every .msg/.srv/.action file to be processed
##   * uncomment the generate_messages entry below
##   * add every package in MSG_DEP_SET to generate_messages(DEPENDENCIES ...)

## Generate messages in the 'msg' folder
# add_message_files(
#   FILES
#   Message1.msg
#   Message2.msg
# )

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
# generate_messages(
#   DEPENDENCIES
#   std_msgs  # Or other packages containing msgs
# )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include
#  LIBRARIES graph_slam
  CATKIN_DEPENDS dynamic_slam_utils ndt_scanmatching2d
  DEPENDS eigen boost pcl
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
# include_directories(include)

## Declare a cpp library
# add_library(graph_slam
#   src/${PROJECT_NAME}/graph_slam.cpp
# )

## Declare a cpp executable
# add_executable(graph_slam_node src/graph_slam_node.cpp)

## Add cmake target dependencies of the executable/library
## as an example, message headers may need to be generated before nodes
# add_dependencies(graph_slam_node graph_slam_generate_messages_cpp)

## Specify libraries to link a library or executable target against
# target_link_libraries(graph_slam_node
#   ${catkin_LIBRARIES}
# )

# c++11 support required
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_C_COMPILER "clang")
set(cxx_compiler_flags "-O3 -msse3 -funroll-loops")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    set(cxx_compiler_flags "${cxx_compiler_flags} -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if(OPENMP_FOUND)
  set(cxx_compiler_flags "${cxx_compiler_flags} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_compiler_flags}")

include_directories(
    include
    ${PCL_INCLUDE_DIRS}
    ${catkin_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${G2O_INCLUDE_DIR}
    ${CSPARSE_INCLUDE_DIR}
)
link_directories(
    ${PCL_LIBRARY_DIRS}
    ${G2O_LIBRARY_DIR}
    ${CSPARSE_LIBRARY_DIRS}
)
add_definitions(
    ${PCL_DEFINITIONS}
)
link_libraries(${G2O_LIBRARIES})

add_executable(graph_slam2d_g2o_parser src/parser2d.cpp)
target_link_libraries(graph_slam2d_g2o_parser
    ${catkin_LIBRARIES}
    ${PCL_COMMON_LIBRARIES}
)

add_executable(graph_slam2d
  src/graph_slam2d_node.cpp
  src/rrr_clusterizer.cpp
  src/rrr_g20_wrapper.cpp
  )
# add_dependencies(graph_slam2d ${${PROJECT_NAME}_EXPORTED_TARGETS})
target_link_libraries(graph_slam2d
   ${catkin_LIBRARIES}
   ${PCL_LIBRARIES}
   ${Boost_LIBRARIES}
   ${G2O_LIBRARIES}
   ${CSPARSE_LIBRARY}
)

add_executable(gauss_newton_optimalizer_test test/gauss_newton_optimalizer_test.cpp)
target_link_libraries(gauss_newton_optimalizer_test
   ${catkin_LIBRARIES}
   ${PCL_LIBRARIES}
   ${Boost_LIBRARIES}
)
#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS graph_slam graph_slam_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

if(CATKIN_ENABLE_TESTING)
  catkin_add_gtest(${PROJECT_NAME}-graph-test
    test/pose_graph_test.cpp
    )
endif()
## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_graph_slam.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
