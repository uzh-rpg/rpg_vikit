/*
 * math_utils.h
 *
 *  Created on: Jul 20, 2012
 *      Author: cforster
 */

#ifndef VIKIT_MATH_UTILS_H_
#define VIKIT_MATH_UTILS_H_

#include <cstdlib> // size_t, fabs
#include <cmath>   // sin, cos
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <sophus/se3.h>
#include <sophus/so3.h>

namespace vk
{

using namespace Eigen;
using Sophus::SE3;
using std::size_t;
using std::uint8_t;

Vector3d triangulateFeatureNonLin(
    const Matrix3d& R,
    const Vector3d& t,
    const Vector3d& feature1,
    const Vector3d& feature2);

/// Assumes the bearing vectors f_c and f_r are on the epipolar plane, i.e.
/// perfect triangulation without noise!
bool depthFromTriangulationExact(
    const Matrix3d& R_r_c,
    const Vector3d& t_r_c,
    const Vector3d& f_r,
    const Vector3d& f_c,
    double& depth_in_r,
    double& depth_in_c);

double reprojError(
    const Vector3d& f1,
    const Vector3d& f2,
    double error_multiplier2);

double computeInliers(
    const std::vector<Vector3d, aligned_allocator<Vector3d> >& features1, ///< c1
    const std::vector<Vector3d, aligned_allocator<Vector3d> >& features2, ///< c2
    const Matrix3d& R, ///< R_c1_c2
    const Vector3d& t, ///< c1_t
    const double reproj_thresh,
    double error_multiplier2,
    std::vector<Vector3d, aligned_allocator<Vector3d> >& xyz_vec, ///< in frame c1
    std::vector<int>& inliers,
    std::vector<int>& outliers
);

void computeInliersOneView(
    const std::vector<Vector3d, aligned_allocator<Vector3d> >& feature_sphere_vec,
    const std::vector<Vector3d, aligned_allocator<Vector3d> >& xyz_vec,
    const Matrix3d &R,
    const Vector3d &t,
    const double reproj_thresh,
    const double error_multiplier2,
    std::vector<int>& inliers,
    std::vector<int>& outliers);

/// Direct Cosine Matrix to Roll Pitch Yaw
Vector3d dcm2rpy(const Matrix3d &R);

/// Roll Pitch Yaw to Direct Cosine Matrix
Matrix3d rpy2dcm(const Vector3d &rpy);

/// Angle Axis parametrization to Quaternion
Quaterniond angax2quat(const Vector3d& n, const double& angle);

/// Angle Axis parametrization to Matrix representation
Matrix3d angax2dcm(const Vector3d& n, const double& angle);

/// Spherical linear interpolation. t should be in [0,1]
Sophus::SO3 slerp(const Sophus::SO3& R0, const Sophus::SO3& R1, double t);

double sampsonusError(
    const Vector2d &v2Dash,
    const Matrix3d& m3Essential,
    const Vector2d& v2);

inline Matrix3d skew(const Vector3d& v)
{
  Matrix3d v_sqew;
  v_sqew << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;
  return v_sqew;
}

/// From GTSAM
inline Matrix3d rightJacobianExpMapSO3(const Vector3d& x)
{
  // x is the axis-angle representation (exponential coordinates) for a rotation
  const double normx = x.norm(); // rotation angle
  Matrix3d Jr;
  if (normx < 10e-8){
    Jr = Matrix3d::Identity();
  }
  else{
    const Matrix3d X = vk::skew(x); // element of Lie algebra so(3): X = x^
    Jr = Matrix3d::Identity() - ((1-cos(normx))/(normx*normx)) * X +
        ((normx-sin(normx))/(normx*normx*normx)) * X * X; // right Jacobian
  }
  return Jr;
}

/// From GTSAM
inline Matrix3d rightJacobianExpMapSO3inverse(const Vector3d& x)
{
  // x is the axis-angle representation (exponential coordinates) for a rotation
  const double normx = x.norm(); // rotation angle
  Matrix3d Jrinv;
  if (normx < 10e-8)
  {
    Jrinv = Matrix3d::Identity();
  }
  else
  {
    const Matrix3d X = vk::skew(x); // element of Lie algebra so(3): X = x^
    Jrinv = Matrix3d::Identity() +
        0.5 * X + (1/(normx*normx) - (1+cos(normx))/(2*normx * sin(normx))   ) * X * X;
  }
  return Jrinv;
}

inline double norm_max(const Eigen::VectorXd & v)
{
  double max = -1;
  for (int i=0; i<v.size(); i++)
  {
    double abs = std::fabs(v[i]);
    if(abs>max){
      max = abs;
    }
  }
  return max;
}

inline Vector2d project2d(const Vector3d& v)
{
  return v.head<2>()/v[2];
}

inline Vector3d unproject2d(const Vector2d& v)
{
  return Vector3d(v[0], v[1], 1.0);
}

inline Vector3d project3d(const Vector4d& v)
{
  return v.head<3>()/v[3];
}

inline Vector4d unproject3d(const Vector3d& v)
{
  return Vector4d(v[0], v[1], v[2], 1.0);
}

template<class T>
T getMedian(std::vector<T>& data_vec)
{
  assert(!data_vec.empty());
  typename std::vector<T>::iterator it = data_vec.begin()+std::floor(data_vec.size()/2);
  std::nth_element(data_vec.begin(), it, data_vec.end());
  return *it;
}

inline double pyrFromZero_d(double x_0, int level)
{
  return x_0/(1<<level); // = 1 / 2^level
}

inline Vector2d pyrFromZero_2d(const Vector2d& uv_0, int level)
{
  return Vector2d(pyrFromZero_d(uv_0[0], level),
                  pyrFromZero_d(uv_0[1], level));
}

/// Frame jacobian for projection of 3D point in (f)rame coordinate to
/// unit plane coordinates uv (focal length = 1).
inline void jacobianFrame_xyz2uv(
    const Vector3d& xyz_in_f,
    Matrix<double,2,6>& J)
{
  const double x = xyz_in_f[0];
  const double y = xyz_in_f[1];
  const double z_inv = 1./xyz_in_f[2];
  const double z_inv_2 = z_inv*z_inv;
  J(0,0) = -z_inv;              // -1/z
  J(0,1) = 0.0;                 // 0
  J(0,2) = x*z_inv_2;           // x/z^2
  J(0,3) = y*J(0,2);            // x*y/z^2
  J(0,4) = -(1.0 + x*J(0,2));   // -(1.0 + x^2/z^2)
  J(0,5) = y*z_inv;             // y/z
  J(1,0) = 0.0;                 // 0
  J(1,1) = -z_inv;              // -1/z
  J(1,2) = y*z_inv_2;           // y/z^2
  J(1,3) = 1.0 + y*J(1,2);      // 1.0 + y^2/z^2
  J(1,4) = -J(0,3);             // -x*y/z^2
  J(1,5) = -x*z_inv;            // x/z
}

/// Jacobian of point projection on unit plane (focal length = 1) in frame (f).
inline void jacobianPoint_xyz2uv(
    const Vector3d& p_in_f,
    const Matrix3d& R_f_w,
    Matrix<double,2,3>& J)
{
  const double z_inv = 1.0/p_in_f[2];
  const double z_inv_sq = z_inv*z_inv;
  J(0,0) = z_inv;
  J(0,1) = 0.0;
  J(0,2) = -p_in_f[0] * z_inv_sq;
  J(1,0) = 0.0;
  J(1,1) = z_inv;
  J(1,2) = -p_in_f[1] * z_inv_sq;
  J = - J * R_f_w;
}

} // end namespace vk

#endif // VIKIT_MATH_UTILS_H_
