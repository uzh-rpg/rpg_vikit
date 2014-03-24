/*
 * blender_utils.h
 *
 *  Created on: Feb 13, 2014
 *      Author: cforster
 */

#ifndef VIKIT_BLENDER_UTILS_H_
#define VIKIT_BLENDER_UTILS_H_

#include <string>
#include <vikit/pinhole_camera.h>
#include <vikit/math_utils.h>
#include <opencv2/opencv.hpp>
#include <fstream>

namespace vk {
namespace blender_utils {

void loadBlenderDepthmap(const std::string file_name,
                         const vk::AbstractCamera& cam,
                         cv::Mat& img)
{
  std::ifstream file_stream(file_name.c_str());
  assert(file_stream.is_open());
  img = cv::Mat(cam.height(), cam.width(), CV_32FC1);
  float * img_ptr = img.ptr<float>();
  float depth;
  for(int y=0; y<cam.height(); ++y)
  {
    for(int x=0; x<cam.width(); ++x, ++img_ptr)
    {
      file_stream >> depth;
      // blender:
      Eigen::Vector2d uv(vk::project2d(cam.cam2world(x,y)));
      *img_ptr = depth * sqrt(uv[0]*uv[0] + uv[1]*uv[1] + 1.0);

      // povray
      // *img_ptr = depth/100.0; // depth is in [cm], we want [m]

      if(file_stream.peek() == '\n' && x != cam.width()-1 && y != cam.height()-1)
        printf("WARNING: did not read the full depthmap!\n");
    }
  }
}

namespace file_format
{

class ImageNameAndPose
{
public:
  ImageNameAndPose() {}
  virtual ~ImageNameAndPose() {}
  double timestamp_;
  std::string image_name_;
  Eigen::Vector3d t_;
  Eigen::Quaterniond q_;
  friend std::ostream& operator <<(std::ostream& out, const ImageNameAndPose& pair);
  friend std::istream& operator >>(std::istream& in, ImageNameAndPose& pair);
};

std::ostream& operator <<(std::ostream& out, const ImageNameAndPose& gt)
{
  out << gt.timestamp_ << " " << gt.image_name_ << " "
      << gt.t_.x()   << " " << gt.t_.y()   << " " << gt.t_.z()   << " "
      << gt.q_.x()   << " " << gt.q_.y()   << " " << gt.q_.z()   << " " << gt.q_.w()   << " " << std::endl;
  return out;
}

std::istream& operator >>(std::istream& in, ImageNameAndPose& gt)
{
  in >> gt.timestamp_;
  in >> gt.image_name_;
  double tx, ty, tz, qx, qy, qz, qw;
  in >> tx;
  in >> ty;
  in >> tz;
  in >> qx;
  in >> qy;
  in >> qz;
  in >> qw;
  gt.t_ = Eigen::Vector3d(tx, ty, tz);
  gt.q_ = Eigen::Quaterniond(qw, qx, qy, qz);
  return in;
}

} // namespace file_format
} // namespace blender_utils
} // namespace vk

#endif // VIKIT_BLENDER_UTILS_H_
