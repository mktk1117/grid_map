/*
 * Registration.cpp
 *
 *  Created on: Dec 4, 2017
 *     Authors: Takahiro Miki
 *   Institute: ispace.inc
 */

#include <grid_map_sdf/SignedDistanceField.hpp>
#include <grid_map_registration/Registration.hpp>
#include <grid_map_pcl/GridMapPclConverter.hpp>

#include <grid_map_core/GridMap.hpp>

 #include <pcl/common/transforms.h>

#include <limits>

using namespace pcl;

namespace grid_map {

Registration::Registration(int maxIteration, double gradientThreshold, double momentThreshold, double distanceThreshold)
  :
    baseMapSet_(false),
    newMapSet_(true)
{
  maxIteration_ = maxIteration;
  gradientThreshold_ = gradientThreshold;
  momentThreshold_ = momentThreshold;
  distanceThreshold_ = distanceThreshold;
}

Registration::~Registration()
{
}


void Registration::setBaseMap(const GridMap& baseMap)
{
  baseMap_ = baseMap;
  baseSdf_.calculateSignedDistanceField(baseMap_, "elevation", 0.5);
  baseSdf_.calculateHorizontalSignedDistanceField(baseMap_, "elevation", 0.5);
  baseMapSet_ = true;
}
void Registration::setNewMap(const GridMap& newMap)
{
  newMap_ = newMap;
  newMapSet_ = true;
}
void Registration::align(GridMap& alignedMap)
{
  Eigen::Affine3d transform;
  align(alignedMap, transform);
}

void Registration::align(GridMap& alignedMap, Eigen::Affine3d& transform)
{
  double resolution = baseMap_.getResolution();
  Position newMapPosition = newMap_.getPosition();
  PointCloud<PointXYZ>::Ptr cloud(new PointCloud<PointXYZ>);
  GridMapPclConverter::toPointCloud(newMap_, "elevation", cloud);

  double cost_th = 10;
  double threshold = 1e-10;

  std::default_random_engine randomGenerator;
  std::normal_distribution<double> randomDistribution(0, resolution * 2);
  Eigen::Vector3d centerOfCloud(0, 0, 0);
  Eigen::Vector3d prevCenterOfCloud(0, 0, 0);

  Eigen::Vector3d translation(0, 0, 0);
  double yaw = 0;
  bool translationConverged = false;

  for (size_t i = 0; i < maxIteration_; i++){
    centerOfCloud << 0, 0, 0;
    Eigen::Vector3d totalGradient(0, 0, 0);
    double totalSquareDistance = 0;
    double totalMoment = 0;
    for (size_t i = 0; i < cloud->points.size (); ++i){
      auto& point = cloud->points[i];
      Eigen::Vector3d p(point.x, point.y, point.z);

      double d = baseSdf_.getInterpolatedDistanceAt(p);
      centerOfCloud += p;
      totalSquareDistance += d * d;
      // signD = 1 if d > 0, -1 if d < 0
      int signD = (d > 0) - (d < 0);

      Eigen::Vector3d gradient = -signD * baseSdf_.getHorizontalDistanceGradientAt(p);
      totalGradient += gradient * fabs(d);
      
      //moment
      if (fabs(d) > resolution && translationConverged){
        Eigen::Vector3d r3d = p - prevCenterOfCloud;
        Eigen::Vector2d r(r3d.x(), r3d.y());
        Eigen::Vector2d g(gradient.x(), gradient.y());
        double moment = (r.x() * g.y() - g.x() * r.y());
        totalMoment += moment * fabs(d);
      }
    }
    totalGradient /= cloud->points.size();
    totalMoment /= cloud->points.size();
    centerOfCloud /= cloud->points.size();
    prevCenterOfCloud = centerOfCloud;

    if (totalSquareDistance < distanceThreshold_ && totalMoment < momentThreshold_)
      break;

    if (totalGradient.norm() < gradientThreshold_){
      translationConverged = true;
      // to avoid local minimum, add random gradient
      if (totalSquareDistance > distanceThreshold_ * 2){
        Eigen::Vector3d r(randomDistribution(randomGenerator), 
            randomDistribution(randomGenerator), randomDistribution(randomGenerator));
        totalGradient += r;
      }
    }

    applyYawRotation(cloud, centerOfCloud, totalMoment);
    applyTranslation(cloud, totalGradient);

    translation += totalGradient;
    yaw += totalMoment;
  }
  transform = Eigen::Affine3d(Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()));
  transform *= Eigen::Translation3d(translation);

  // make aligned map
  alignedMap = baseMap_;
  GridMapPclConverter::toGridMap(alignedMap, "elevation", cloud);
}

void Registration::applyYawRotation(PointCloud<PointXYZ>::Ptr cloud, Eigen::Vector3d center, double yaw)
{
  Eigen::Affine3d transform;
  transform = Eigen::Affine3d(Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()));
  for (size_t i = 0; i < cloud->points.size (); ++i){
    cloud->points[i].x -= center.x();
    cloud->points[i].y -= center.y();
    cloud->points[i].z -= center.z();
  }
  pcl::transformPointCloud(*cloud, *cloud, transform.cast<float>());
  for (size_t i = 0; i < cloud->points.size (); ++i){
    cloud->points[i].x += center.x();
    cloud->points[i].y += center.y();
    cloud->points[i].z += center.z();
  }
}

void Registration::applyTranslation(PointCloud<PointXYZ>::Ptr cloud, Eigen::Vector3d translation)
{
  for (size_t i = 0; i < cloud->points.size (); ++i){
    cloud->points[i].x += translation.x();
    cloud->points[i].y += translation.y();
    cloud->points[i].z += translation.z();
  }
}




} /* namespace */
