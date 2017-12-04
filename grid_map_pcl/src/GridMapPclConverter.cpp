/*
 * GridMapPclConverter.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: Dominic Jud
 *   Institute: ETH Zurich, Robotic Systems Lab
 */

#include "grid_map_pcl/GridMapPclConverter.hpp"

namespace grid_map {

GridMapPclConverter::GridMapPclConverter()
{
}

GridMapPclConverter::~GridMapPclConverter()
{
}

bool GridMapPclConverter::initializeFromPolygonMesh(const pcl::PolygonMesh& mesh,
                                                    const double resolution,
                                                    grid_map::GridMap& gridMap)
{
  pcl::PointCloud < pcl::PointXYZ > cloud;
  pcl::fromPCLPointCloud2(mesh.cloud, cloud);
  pcl::PointXYZ minBound;
  pcl::PointXYZ maxBound;
  pcl::getMinMax3D(cloud, minBound, maxBound);

  grid_map::Length length = grid_map::Length(maxBound.x - minBound.x, maxBound.y - minBound.y);
  grid_map::Position position = grid_map::Position((maxBound.x + minBound.x) / 2.0,
                                                   (maxBound.y + minBound.y) / 2.0);
  gridMap.setGeometry(length, resolution, position);

  return true;
}

bool GridMapPclConverter::addLayerFromPolygonMesh(const pcl::PolygonMesh& mesh,
                                                  const std::string& layer,
                                                  grid_map::GridMap& gridMap)
{
  // Adding a layer to the grid map to put data into
  gridMap.add(layer);
  // Converting out of binary cloud data
  pcl::PointCloud <pcl::PointXYZ> cloud;
  pcl::fromPCLPointCloud2(mesh.cloud, cloud);
  // Direction and max height for projection ray
  const Eigen::Vector3f ray = -Eigen::Vector3f::UnitZ();
  pcl::PointXYZ minBound;
  pcl::PointXYZ maxBound;
  pcl::getMinMax3D(cloud, minBound, maxBound);

  // Iterating over the triangles in the mesh
  for (const pcl::Vertices& polygon : mesh.polygons) {
    // Testing this is a triangle
    assert(polygon.vertices.size() == 3);
    // Getting the vertices of the triangle (as a single matrix)
    Eigen::Matrix3f triangleVertexMatrix;
    triangleVertexMatrix.row(0) = cloud[polygon.vertices[0]].getVector3fMap();
    triangleVertexMatrix.row(1) = cloud[polygon.vertices[1]].getVector3fMap();
    triangleVertexMatrix.row(2) = cloud[polygon.vertices[2]].getVector3fMap();
    // Getting the bounds in the XY plane (down projection)
    float maxX = triangleVertexMatrix.col(0).maxCoeff();
    float minX = triangleVertexMatrix.col(0).minCoeff();
    float maxY = triangleVertexMatrix.col(1).maxCoeff();
    float minY = triangleVertexMatrix.col(1).minCoeff();
    // Iterating over the grid cells in the a submap below the triangle
    grid_map::Length length(maxX - minX, maxY - minY);
    grid_map::Position position((maxX + minX) / 2.0, (maxY + minY) / 2.0);
    bool isSuccess;
    SubmapGeometry submap(gridMap, position, length, isSuccess);
    if (isSuccess) {
      for (grid_map::SubmapIterator iterator(submap); !iterator.isPastEnd();
           ++iterator) {
        // Cell position
        const Index index(*iterator);
        grid_map::Position vertexPositionXY;
        gridMap.getPosition(index, vertexPositionXY);
        // Ray origin
        Eigen::Vector3f point(vertexPositionXY.x(), vertexPositionXY.y(),
                              maxBound.z + 1.0);
        // Vertical ray/triangle intersection
        Eigen::Vector3f intersectionPoint;
        if (rayTriangleIntersect(point, ray, triangleVertexMatrix, intersectionPoint)) {
          // If data already present in this cell, taking the max, else setting the data
          if (gridMap.isValid(index, layer)) {
            gridMap.at(layer, index) = std::max(gridMap.at(layer, index), intersectionPoint.z());
          } else {
            gridMap.at(layer, index) = intersectionPoint.z();
          }
        }
      }
    }
  }
  // Success
  return true;
}

bool GridMapPclConverter::rayTriangleIntersect(
    const Eigen::Vector3f& point, const Eigen::Vector3f& ray,
    const Eigen::Matrix3f& triangleVertexMatrix,
    Eigen::Vector3f& intersectionPoint) {
  // Algorithm here is adapted from:
  // http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
  //
  // Original copyright notice:
  // Copyright 2001, softSurfer (www.softsurfer.com)
  // This code may be freely used and modified for any purpose
  // providing that this copyright notice is included with it.

  const Eigen::Vector3f a = triangleVertexMatrix.row(0);
  const Eigen::Vector3f b = triangleVertexMatrix.row(1);
  const Eigen::Vector3f c = triangleVertexMatrix.row(2);
  const Eigen::Vector3f u = b - a;
  const Eigen::Vector3f v = c - a;
  const Eigen::Vector3f n = u.cross(v);
  const float n_dot_ray = n.dot(ray);

  if (std::fabs(n_dot_ray) < 1e-9) return false;

  const float r = n.dot(a - point) / n_dot_ray;

  if (r < 0) return false;

  // Note(alexmillane): The addition of this comparison delta (not in the
  // original algorithm) means that rays intersecting the edge of triangles are
  // treated at hits.
  constexpr float delta = 1e-5;

  const Eigen::Vector3f w = point + r * ray - a;
  const float denominator = u.dot(v) * u.dot(v) - u.dot(u) * v.dot(v);
  const float s_numerator = u.dot(v) * w.dot(v) - v.dot(v) * w.dot(u);
  const float s = s_numerator / denominator;
  if (s < (0 - delta) || s > (1 + delta)) return false;

  const float t_numerator = u.dot(v) * w.dot(u) - u.dot(u) * w.dot(v);
  const float t = t_numerator / denominator;
  if (t < (0 - delta) || s + t > (1 + delta)) return false;

  intersectionPoint = a + s * u + t * v;

  return true;
}

void GridMapPclConverter::toPointCloud(const grid_map::GridMap& gridMap,
                         const std::string& pointLayer,
                         pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud)
{
  double resolution = gridMap.getResolution();
  for (GridMapIterator it(gridMap); !it.isPastEnd(); ++it) {
    Position position;
    if(!gridMap.getPosition(*it, position)) continue;
    double centerValue = gridMap.at(pointLayer, *it);
    double maxValue = std::numeric_limits<float>::min();
    if (std::isnan(centerValue)) continue;
    for (int i = -1; i <= 1; i++){
      for (int j = -1; j <= 1; j++){
        if( i * j == 0 && !(i == 0 && j == 0)){
          Position p = position + Position(i * resolution, j * resolution);
          Eigen::Array2i index_p;
          if(!gridMap.getIndex(p, index_p)) continue;
          double value = gridMap.at(pointLayer, index_p);
          if( value < centerValue && value > maxValue){
              maxValue = value;
          }
        }
      }
    }
    pointCloud->push_back(pcl::PointXYZ(position.x(), position.y(), centerValue));
    if( maxValue > std::numeric_limits<float>::min()){
      double valueDifference = centerValue - maxValue;
      for(int i = 0; i < valueDifference / resolution; i++){
        double z = centerValue - resolution * i;
        pcl::PointXYZ p(position.x(), position.y(), z);
        pointCloud->push_back(p);
      }
    }
  }
}

void GridMapPclConverter::toGridMap(grid_map::GridMap& gridMap,
                         const std::string& pointLayer,
                         const pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud)
{
  for (size_t i = 0; i < pointCloud->points.size (); ++i){
    const auto& point = pointCloud->points[i];
    Position position(point.x, point.y);
    Index index;
    if(!gridMap.getIndex(position, index)) continue;
    auto& elevation = gridMap.at("elevation", index);
    if (std::isnan(elevation) || point.z > elevation)
      elevation = point.z;
  }
}

} /* namespace */
