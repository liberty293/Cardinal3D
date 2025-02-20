
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.
    float xmin, xmax, ymin, ymax, zmin, zmax;
    if (ray.dir.x == 0) {
        if (min.x <= ray.point.x && ray.point.x <= max.x)
            xmin = times.x, xmax = times.y;
        else
            return false;
    } else {
        xmin = (min.x - ray.point.x) / ray.dir.x;
        xmax = (max.x - ray.point.x) / ray.dir.x;
        if (xmin > xmax)
          std::swap(xmin, xmax);
    }

    if (ray.dir.y == 0) {
        if (min.y <= ray.point.y && ray.point.y <= max.y)
            ymin = times.x, ymax = times.y;
        else
            return false;
    } else {
        ymin = (min.y - ray.point.y) / ray.dir.y;
        ymax = (max.y - ray.point.y) / ray.dir.y;
        if (ymin > ymax)
          std::swap(ymin, ymax);
    }

    if (ray.dir.z == 0) {
        if (min.z <= ray.point.z && ray.point.z <= max.z)
            zmin = times.x, zmax = times.y;
        else
            return false;
    } else {
        zmin = (min.z - ray.point.z) / ray.dir.z;
        zmax = (max.z - ray.point.z) / ray.dir.z;
        if (zmin > zmax)
          std::swap(zmin, zmax);
    }

    float amin = std::max(std::max(std::max(ray.dist_bounds.x, times.x), xmin), std::max(ymin, zmin));
    float amax = std::min(std::min(std::min(ray.dist_bounds.y, times.y), xmax), std::min(ymax, zmax));

    times.x = amin, times.y = amax;
    return amin <= amax;
}
