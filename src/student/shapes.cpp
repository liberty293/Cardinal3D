
#include "../rays/shapes.h"
#include "debug.h"

namespace PT {

const char* Shape_Type_Names[(int)Shape_Type::count] = {"None", "Sphere"};

BBox Sphere::bbox() const {

    BBox box;
    box.enclose(Vec3(-radius));
    box.enclose(Vec3(radius));
    return box;
}

Trace Sphere::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!
    Trace ret;
    ret.origin = ray.point;
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?

    //implcid surface all |x|^2 = radius ^2
    float b = dot(ray.point, ray.dir);
    float c = ray.point.norm_squared()-(radius*radius);
    if(b*b-c < 0)
        return ret; //no intersection

    float t1 = -b + sqrt(b*b-c);
    float t2 = -b - sqrt(b*b-c);
    float tmin = t1 < t2 ? t1:t2;
    float tmax = tmin==t1? t2:t1;
    if(tmin < ray.dist_bounds.y && tmin > ray.dist_bounds.x) //closest point and in bounds
    {
        ret.hit = true;
        ret.distance = tmin;
        ret.position = ray.point + ray.dir*tmin;
        ret.normal = ret.position.normalize(); //from the origin to the point out of sphere
        ray.dist_bounds.y = tmin;
        return ret;
    }
    else if(tmax < ray.dist_bounds.y && tmax > ray.dist_bounds.x) //if close isn't in bounds check far
    {
        ret.hit = true;
        ret.distance = tmax;
        ret.position = ray.point + ray.dir*tmax;
        ret.normal = ret.position.normalize(); //from the origin to the point out of sphere
        ray.dist_bounds.y = tmax;
        return ret;
    }
    

    return ret;
}

} // namespace PT
