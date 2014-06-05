cmake_minimum_required(VERSION 2.8)
include($ENV{LUCY_ROOT}/cmake/LucyCommon.cmake)

setup_module()

##
## VIKIT_COMMON
##

set(SOURCES
   vikit_common/src/img_align.cpp
   vikit_common/src/homography.cpp
   vikit_common/src/pinhole_camera.cpp
   vikit_common/src/robust_cost.cpp
   vikit_common/src/vision.cpp
   vikit_common/src/user_input_thread.cpp
   vikit_common/src/omni_camera.cpp
   vikit_common/src/performance_monitor.cpp
   vikit_common/src/math_utils.cpp
   vikit_common/src/atan_camera.cpp
   )

set(HEADERS
   vikit_common/include/vikit/omni_camera.h
   vikit_common/include/vikit/user_input_thread.h
   vikit_common/include/vikit/sample.h
   vikit_common/include/vikit/robust_cost.h
   vikit_common/include/vikit/performance_monitor.h
   vikit_common/include/vikit/file_reader.h
   vikit_common/include/vikit/aligned_mem.h
   vikit_common/include/vikit/img_align.h
   vikit_common/include/vikit/patch_score.h
   vikit_common/include/vikit/vision.h
   vikit_common/include/vikit/blender_utils.h
   vikit_common/include/vikit/nlls_solver_impl.hpp
   vikit_common/include/vikit/atan_camera.h
   vikit_common/include/vikit/ringbuffer.h
   vikit_common/include/vikit/nlls_solver.h
   vikit_common/include/vikit/pinhole_camera.h
   vikit_common/include/vikit/abstract_camera.h
   vikit_common/include/vikit/timer.h
   vikit_common/include/vikit/math_utils.h
   vikit_common/include/vikit/homography.h
   )

message(STATUS "src_dir=${PROJECT_SOURCE_DIR}/vikit_common/include")

create_module()

add_to_system_path(${PROJECT_SOURCE_DIR}/vikit_common/include)

include_3rdparty(sophus)

find_boost()
find_eigen3()
