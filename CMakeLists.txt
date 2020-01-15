cmake_minimum_required(VERSION 3.14)
project(PointLinePlaneParam)


SET(CMAKE_BUILD_TYPE Release) # Release, RelWithDebInfo , Debug
set(CMAKE_CXX_STANDARD 14)

find_package(Pangolin REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

INCLUDE_DIRECTORIES(
        ${CMAKE_SOURCE_DIR}
        ${EIGEN3_INCLUDE_DIR}
        ${Pangolin_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS}
        ${CERES_INCLUDE_DIRS}
        include
)

file(GLOB SOURCEFILES src/*.cc include/*.hpp)

add_library(PointLinePlaneParamLib SHARED ${SOURCEFILES})

target_link_libraries(PointLinePlaneParamLib
        ${OpenCV_LIBS}
        ${EIGEN3_LIBS}
        ${Pangolin_LIBRARIES}
        ${CERES_LIBRARIES}
        )

add_executable(test_pointplane test/test_pointplane.cc)
target_link_libraries(test_pointplane PointLinePlaneParamLib ${CMAKE_THREAD_LIBS_INIT})