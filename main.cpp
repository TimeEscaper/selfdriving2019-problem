/**
 * Solution for problem https://contest.yandex.ru/contest/12698/problems/
 *
 * Solution is based on RANSAC approach for plane fitting in point cloud.
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>

/**
 * Point descriptor.
 */
typedef struct
{
  double x;
  double y;
  double z;
} Point3d;

/**
 * Plane descriptor.
 */
typedef struct
{
  double a;
  double b;
  double c;
  double d;
} Plane3d;

/**
 * Creates plane from 3 points.
 *
 * @param point1 First point
 * @param point2 Second point
 * @param point3 Third point
 *
 * @return Plane descriptor
 */
Plane3d createPlane(const Point3d& point1, const Point3d& point2, const Point3d& point3);

/**
 * Picks a random sub-vector from specified vector.
 *
 * @param vector Vector to pick from
 * @param count Number of elements to pick
 *
 * @return Vector of size {@param count} with random elements from original {@param vector}
 */
std::vector<Point3d> getRandomSubVector(const std::vector<Point3d>& vector, long count);

/**
 * Calculated distance from point to plane.
 *
 * @param point Point for calculations
 * @param plane Plane for calculations
 *
 * @return Distance from {@param point} to {@param plane}
 */
double getDistanceToPlane(const Point3d& point, const Plane3d& plane);

/**
 * Input file name.
 */
static const std::string INPUT_FILE = "input.txt";

/**
 * Maximal number of iterations for RANSAC approach.
 */
static const int RANSAC_MAX_ITERATIONS = 100;

/**
 * Number of elements for one RANSAC hypothesis.
 */
static const int RANSAC_NUM_ELEMENTS = 3;

int main()
{
  double p = 0.0;
  long N = 0;
  std::vector<Point3d> points;

  // Read data from file
  std::ifstream inputStream(INPUT_FILE.c_str());
  inputStream >> p;
  inputStream >> N;
  double xBuf, yBuf, zBuf = 0.0;
  Point3d pointBuf;
  while (inputStream >> xBuf >> yBuf >> zBuf)
  {
    pointBuf.x = xBuf; pointBuf.y = yBuf; pointBuf.z = zBuf;
    points.push_back(pointBuf);
  }

  // RANSAC parameters
  int iterations = RANSAC_MAX_ITERATIONS;
  int inlinersThreshold = points.size() / 2.0;

  Plane3d bestPlane;
  long bestInliners = 0;

  while (iterations--)
  {
    std::vector<Point3d> samples = getRandomSubVector(points, RANSAC_NUM_ELEMENTS);
    Plane3d plane = createPlane(samples[0], samples[1], samples[2]);

    long inliners = 0;
    for (auto & point : points)
    {
      if (getDistanceToPlane(point, plane) <= p)
      {
        inliners++;
      }
    }

    if (inliners > bestInliners)
    {
      bestInliners = inliners;
      bestPlane = plane;
    }

    if (inliners > inlinersThreshold)
    {
      break;
    }
  }

  std::cout << bestPlane.a << " " << bestPlane.b << " " << bestPlane.c << " " << bestPlane.d << std::endl;

  return 0;
}

std::vector<Point3d> getRandomSubVector(const std::vector<Point3d>& vector, long count)
{
  std::set<long> chosen;
  std::vector<Point3d> result;
  while (count != 0)
  {
    auto index = std::rand() % vector.size();
    if (chosen.find(index) != chosen.end())
    {
      continue;
    }
    result.push_back(vector[index]);
    chosen.insert(index);
    count--;
  }
  return result;
}

double getDistanceToPlane(const Point3d& point, const Plane3d& plane)
{
  return std::abs(plane.a*point.x + plane.b*point.y + plane.c*point.z + plane.d) /
    std::sqrt(plane.a*plane.a + plane.b*plane.b + plane.c*plane.c);
}

Plane3d createPlane(const Point3d& point1, const Point3d& point2, const Point3d& point3)
{
  Plane3d result;
  double a1 = point2.x - point1.x;
  double b1 = point2.y - point1.y;
  double c1 = point2.z - point1.z;
  double a2 = point3.x - point1.x;
  double b2 = point3.y - point1.y;
  double c2 = point3.z - point1.z;
  double a = b1 * c2 - b2 * c1;
  double b = a2 * c1 - a1 * c2;
  double c = a1 * b2 - b1 * a2;
  double d = (- a * point1.x - b * point1.y - c * point1.z);
  double norm = std::sqrt(a*a + b*b + c*c + d*d);
  result.a = a / norm; result.b = b / norm; result.c = c / norm; result.d = d / norm;
  return result;
}
