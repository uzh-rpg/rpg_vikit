/*
 * homography.cpp
 * Adaptation of PTAM-GPL HomographyInit class.
 * https://github.com/Oxford-PTAM/PTAM-GPL
 * Licence: GPLv3
 * Copyright 2008 Isis Innovation Limited
 *
 *  Created on: Sep 2, 2012
 *      by: cforster
 *
 * This class implements the homography decomposition of Faugeras and Lustman's
 * 1988 tech report. Code converted to Eigen from PTAM.
 *
 */

#ifndef VIKIT_HOMOGRAPHY_H_
#define VIKIT_HOMOGRAPHY_H_

#include <vikit/math_utils.h>
#include <Eigen/SVD>

namespace vk {

struct HomographyDecomposition
{
  Vector3d t;
  Matrix3d R;
  double d;
  Vector3d n;

  // Resolved  Composition
  Sophus::SE3 T; //!< second from first
  int score;
};

class Homography
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Homography(
         const std::vector<Vector2d, aligned_allocator<Vector2d> >& _fts1,
         const std::vector<Vector2d, aligned_allocator<Vector2d> >& _fts2,
         double error_multiplier2,
         double thresh_in_px);

  void calcFromPlaneParams(
          const Vector3d& normal,
          const Vector3d& point_on_plane);

  void calcFromMatches();

  size_t computeMatchesInliers();

  bool computeSE3fromMatches();

  bool decompose();

  void findBestDecomposition ();

  double thresh;
  double error_multiplier2;
  const std::vector<Vector2d, aligned_allocator<Vector2d> >& fts_c1; //!< Features on first image on unit plane
  const std::vector<Vector2d, aligned_allocator<Vector2d> >& fts_c2; //!< Features on second image on unit plane
  std::vector<bool> inliers;
  SE3 T_c2_from_c1;             //!< Relative translation and rotation of two images
  Matrix3d H_c2_from_c1;                   //!< Homography
  std::vector<HomographyDecomposition> decompositions;
};

} // namespace vk

#endif // VIKIT_HOMOGRAPHY_H_
