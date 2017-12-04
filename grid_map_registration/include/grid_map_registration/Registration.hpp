/*
 * Registration.hpp
 *
 *  Created on: Dec 4, 2017
 *     Authors: Takahiro Miki
 *   Institute: ispace.inc
 */

#pragma once

#include <grid_map_core/GridMap.hpp>
#include <grid_map_sdf/SignedDistanceField.hpp>

#include <pcl/point_types.h>
#include <pcl/conversions.h>

#include <string>
#include <vector>

using namespace pcl;

namespace grid_map {

class Registration
{
 public:
  Registration(int maxIteration, double gradientThreshold, double momentThreshold, double distanceThreshold);
  virtual ~Registration();

  void setBaseMap(const GridMap& baseMap);
  void setNewMap(const GridMap& newMap);
  void align(GridMap& alignedMap);
  void align(GridMap& alignedMap, Eigen::Affine3d& transform);



 private:
  GridMap baseMap_;
  GridMap newMap_;
  bool baseMapSet_;
  bool newMapSet_;

  int maxIteration_;
  double gradientThreshold_;
  double momentThreshold_;
  double distanceThreshold_;

  SignedDistanceField baseSdf_;

  void applyYawRotation(PointCloud<PointXYZ>::Ptr cloud, Eigen::Vector3d center, double yaw);
  void applyTranslation(PointCloud<PointXYZ>::Ptr cloud, Eigen::Vector3d translation);
};

} /* namespace */
