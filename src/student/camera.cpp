
#include "../util/camera.h"
#include "../rays/samplers.h"
#include "debug.h"

Ray Camera::generate_ray(Vec2 screen_coord) const {

    // TODO (PathTracer): Task 1
    //
    // The input screen_coord is a normalized screen coordinate [0,1]^2
    //
    // You need to transform this 2D point into a 3D position on the sensor plane, which is
    // located one unit away from the pinhole in camera space (aka view space).
    Vec3 pos(0,0,-1.0);
    //
    // You'll need to compute this position based on the vertial field of view
    // (vert_fov) of the camera, and the aspect ratio of the output image (aspect_ratio).
    float vheight = tan(vert_fov*(PI_F/180)/2)*2*focal_dist; //1 unit away
    float width = vheight*aspect_ratio;
    Vec2 ratio(width, vheight);
    Vec2 xy = screen_coord*ratio-ratio/2;
    pos.x += xy.x;
    pos.y += xy.y;
    //
    // Tip: compute the ray direction in view space in homogeneous w=0) and use
    // the camera space to world space transform (iview) to transform the ray back into world space.
    Vec3 o(iview*Vec3(0));

    return Ray(o, iview*pos-o);
}
