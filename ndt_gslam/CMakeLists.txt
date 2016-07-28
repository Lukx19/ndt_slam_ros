cmake_minimum_required(VERSION 2.8.3)
project(ndt_gslam)

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
    tf_conversions
    pcl_ros
    message_generation
)

find_package(PCL 1.7 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options system)
find_package(Eigen3 REQUIRED)
# g2o #######################
find_package(G2O REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED )
find_package(CSparse REQUIRED )
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
add_message_files(
  FILES
  NDTMapMsg.msg
  NDTCellMsg.msg
)

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
generate_messages(
  DEPENDENCIES
  geometry_msgs
  std_msgs
)

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
  CATKIN_DEPENDS message_runtime
  DEPENDS eigen boost libpcl-all-dev opencv2
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
set(cxx_compiler_flags "-mfpmath=sse -msse3 -fPIC -O3")
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
    ${OpenCV_INCLUDE_DIRS}
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
set(MAX_MIXTURES_SRC_FILES
  src/slam_optimizer/max_mixture/edge_se2_mixture.cpp
  src/slam_optimizer/max_mixture/edge_se2_pointxy_mixture.cpp
  src/slam_optimizer/max_mixture/edge_se3_mixture.cpp
  src/slam_optimizer/max_mixture/types_g2o_mixture.cpp
  )

set(SWITCH_CONSTRAINS_SRC_FILES
  src/slam_optimizer/switchable_constraints/edge_se2Switchable.cpp
  # src/slam_optimizer/switchable_constraints/edge_se3Switchable.cpp
  src/slam_optimizer/switchable_constraints/edge_switchPrior.cpp
  src/slam_optimizer/switchable_constraints/vertex_switchLinear.cpp
  src/slam_optimizer/switchable_constraints/types_g2o_robust.cpp
)

add_executable(graph_slam2d
  src/graph_slam2d_node.cpp
  src/ndt/ndt_cell.cpp
  ${MAX_MIXTURES_SRC_FILES}
  ${SWITCH_CONSTRAINS_SRC_FILES}
  )
add_dependencies(graph_slam2d ${PROJECT_NAME}_generate_messages_cpp)
  #${${PROJECT_NAME}_EXPORTED_TARGETS})
target_link_libraries(graph_slam2d
   ${catkin_LIBRARIES}
   ${PCL_LIBRARIES}
   ${Boost_LIBRARIES}
   ${G2O_LIBRARIES}
   ${CSPARSE_LIBRARY}
   ${OpenCV_LIBS}
   ${EIGEN3_LIBRARIES}
)

add_executable(ndt2d_visual src/registration/ndt2d_visual.cpp)
target_link_libraries (ndt2d_visual
  ${PCL_LIBRARIES}
  ${catkin_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  )
add_executable(ndt_grid2d_merging
  src/ndt/ndt_grid2d_merging.cpp
  src/ndt/ndt_cell.cpp
  ${MAX_MIXTURES_SRC_FILES}
  ${SWITCH_CONSTRAINS_SRC_FILES}
  )
target_link_libraries (ndt_grid2d_merging
  ${PCL_LIBRARIES}
  ${catkin_LIBRARIES}
  ${OpenCV_LIBS}
  ${EIGEN3_LIBRARIES}
  ${G2O_LIBRARIES}
  ${CSPARSE_LIBRARY}
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
#install(TARGETS ${PROJECT_NAME}
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
#)

## Mark cpp header files for installation
install(DIRECTORY include/${PROJECT_NAME}/
   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
 )

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
    test/slam_optimizer/pose_graph_test.cpp
    )

  catkin_add_gtest(${PROJECT_NAME}-eigen-tools-test
    test/utils/eigen_tools_test.cpp)

  catkin_add_gtest(${PROJECT_NAME}-voxel_grid2d-test
    test/ndt/voxel_grid2d_test.cpp)

  # catkin_add_gtest(${PROJECT_NAME}-ndt2d-test
  #   test/registration/ndt2d_test.cpp
  #   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
  #   )
  if(TARGET ${PROJECT_NAME}-ndt2d-test)
   target_link_libraries(${PROJECT_NAME}-ndt2d-test
    ${catkin_LIBRARIES}
    ${PCL_COMMON_LIBRARIES}
    ${PCL_IO_LIBRARIES}
    ${PCL_REGISTRATION_LIBRARIES}
    ${Boost_LIBRARIES}
    )
  endif()
endif()
## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_graph_slam.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
