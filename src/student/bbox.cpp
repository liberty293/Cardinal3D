
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.
    //check the x coordinate
    if(ray.dir.x != 0) //parallel with x axis; will never hit
    {
        float ax = 1/ray.dir.x;
        float bx = -ray.point.x/ray.dir.x;
        Vec2 tx(ax*min.x +bx, ax*max.x + bx);
        if(times.x > tx.y || times.y  < tx.x) //outside of the times
            return false;
        times.x = times.x > tx.x ? times.x : tx.x; //take the greatest of the small range
        times.y = times.y < tx.y ? times.y : tx.y; //and the least of the big range
    }
    else
    {
        if(ray.point.x > max.x || ray.point.x < min.x) //not in between the planes
            return false;
    }

    //check y coordinate
    if(ray.dir.y != 0) //parallel with y axis; will never hit
    {
        float a = 1/ray.dir.y;
        float b = -ray.point.y/ray.dir.y;
        Vec2 t(a*min.y +b, a*max.y + b);
        if(times.x > t.y || times.y  < t.x) //outside of the times
            return false;
        times.x = times.x > t.x ? times.x : t.x; //take the greatest of the small range
        times.y = times.y < t.y ? times.y : t.y; //and the least of the big range
    }
    else
    {
        if(ray.point.y > max.y || ray.point.y < min.y) //not in between the planes
            return false;
    }
    //check z coordinate
    if(ray.dir.z != 0) //parallel with y axis; will never hit
    {
        float a = 1/ray.dir.z;
        float b = -ray.point.z/ray.dir.z;
        Vec2 t(a*min.z +b, a*max.z + b);
        if(times.x > t.y || times.y  < t.x) //outside of the times
            return false;
        times.x = times.x > t.x ? times.x : t.x; //take the greatest of the small range
        times.y = times.y < t.y ? times.y : t.y; //and the least of the big range
    }
    else
    {
        if(ray.point.z > max.z || ray.point.z < min.z) //not in between the planes
            return false;
    }

    return true;
}
