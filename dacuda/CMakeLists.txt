SET(PROJECT_NAME vikit_common)

if(DEFINED ENV{HOST_LIBS})
   set(CMAKE_INSTALL_PREFIX $ENV{HOST_LIBS}/rpg_vikit CACHE PATH "Install path prefix" FORCE)
endif()

PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.3)
SET(CMAKE_BUILD_TYPE Release) # Release, RelWithDebInfo
SET(CMAKE_VERBOSE_MAKEFILE OFF)
#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

### TODO: enable ARM support if we want to have this for Android...

# # Set IS_ARM on odroid board as environment variable
# IF(DEFINED ENV{ARM_ARCHITECTURE})
#   SET(IS_ARM TRUE)
# ELSE()
SET(IS_ARM FALSE)
# ENDIF()

# Set build flags
set(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS "-Wall -D_LINUX -D_REENTRANT -march=native -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-unused-local-typedefs")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
IF(IS_ARM)
   SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -mfpu=neon -march=armv7 -fomit-frame-pointer -fno-signed-zeros -fno-math-errno -funroll-loops")
ELSE()
   SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -mmmx -msse -msse -msse2 -msse3 -mssse3 -fomit-frame-pointer -fsee -fno-signed-zeros -fno-math-errno -funroll-loops")
ENDIF()
# IF(CMAKE_COMPILER_IS_GNUCC)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
# ELSE()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# ENDIF()

# # Add catkin and required ROS packages
# FIND_PACKAGE(catkin REQUIRED COMPONENTS roscpp cmake_modules)

# Add plain cmake packages 
FIND_PACKAGE(OpenCV REQUIRED)
FIND_PACKAGE(Eigen REQUIRED)
FIND_PACKAGE(Sophus REQUIRED)

# TODO (MWE) handle boost dependent on your setup
add_definitions(-DBOOST_ALL_NO_LIB)
if(LUCY_SYSTEM_WINDOWS)
   set(Boost_USE_STATIC_LIBS ON)
endif()

if(DEFINED ENV{BOOST_ROOT})
   # cmake is too stupid to check for the set environment variable
   set(Boost_NO_SYSTEM_PATHS ON)
elseif(DEFINED ENV{HOST_LIBS})
   #set(Boost_NO_SYSTEM_PATHS ON)
   set(BOOST_ROOT $ENV{HOST_LIBS}/boost)
   message(STATUS "BOOST_ROOT=${BOOST_ROOT}")
endif()
set(Boost_USE_MULTITHREADED ON)

FIND_PACKAGE(Boost REQUIRED COMPONENTS signals system thread)


set(VIKIT_COMMON_RELDIR ../vikit_common)

# Include dirs
INCLUDE_DIRECTORIES(
   ${VIKIT_COMMON_RELDIR}/include
   ${Eigen_INCLUDE_DIRS}
   ${OpenCV_INCLUDE_DIRS}
   ${Sophus_INCLUDE_DIRS}
   ${Boost_INCLUDE_DIRS}
   ${catkin_INCLUDE_DIRS}
   )

# # Describe catkin Project
# catkin_package(
#    DEPENDS Eigen OpenCV Sophus Boost
#    CATKIN_DEPENDS roscpp
#    INCLUDE_DIRS include
#    LIBRARIES ${PROJECT_NAME}
#    )

# Set Sourcefiles
LIST(APPEND SOURCEFILES
   ${VIKIT_COMMON_RELDIR}/src/atan_camera.cpp
   ${VIKIT_COMMON_RELDIR}/src/omni_camera.cpp
   ${VIKIT_COMMON_RELDIR}/src/math_utils.cpp
   ${VIKIT_COMMON_RELDIR}/src/vision.cpp
   ${VIKIT_COMMON_RELDIR}/src/performance_monitor.cpp
   ${VIKIT_COMMON_RELDIR}/src/robust_cost.cpp
   ${VIKIT_COMMON_RELDIR}/src/user_input_thread.cpp
   ${VIKIT_COMMON_RELDIR}/src/pinhole_camera.cpp
   ${VIKIT_COMMON_RELDIR}/src/homography.cpp
   ${VIKIT_COMMON_RELDIR}/src/img_align.cpp
)

LIST(APPEND HEADERS
   ${VIKIT_COMMON_RELDIR}/include/vikit/omni_camera.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/user_input_thread.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/sample.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/robust_cost.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/performance_monitor.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/file_reader.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/aligned_mem.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/img_align.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/patch_score.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/vision.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/blender_utils.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/nlls_solver_impl.hpp
   ${VIKIT_COMMON_RELDIR}/include/vikit/atan_camera.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/ringbuffer.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/nlls_solver.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/pinhole_camera.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/abstract_camera.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/timer.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/math_utils.h
   ${VIKIT_COMMON_RELDIR}/include/vikit/homography.h
   )

# Create vikit library
ADD_LIBRARY(${PROJECT_NAME} STATIC ${SOURCEFILES})
TARGET_LINK_LIBRARIES(
   ${PROJECT_NAME}
   ${OpenCV_LIBS}
   ${Sophus_LIBRARIES}
)#   ${catkin_LIBRARIES})

INSTALL(TARGETS ${PROJECT_NAME}
   #RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
   LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
   ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/vikit)


# # Tests
# ADD_EXECUTABLE(test_vk_common_camera test/test_camera.cpp)
# TARGET_LINK_LIBRARIES(test_vk_common_camera ${PROJECT_NAME})

# ADD_EXECUTABLE(test_vk_common_triangulation test/test_triangulation.cpp)
# TARGET_LINK_LIBRARIES(test_vk_common_triangulation ${PROJECT_NAME})

# ADD_EXECUTABLE(test_vk_common_patch_score test/test_patch_score.cpp)
# TARGET_LINK_LIBRARIES(test_vk_common_patch_score ${PROJECT_NAME})
