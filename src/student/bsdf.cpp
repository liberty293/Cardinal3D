
#include "../rays/bsdf.h"
#include "../util/rand.h"
#include "debug.h"

namespace PT {

Vec3 reflect(Vec3 dir) {

    // TODO (PathTracer): Task 6
    // Return reflection of dir about the surface normal (0,1,0).
    return Vec3(-dir.x, dir.y, -dir.z);
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {

    // TODO (PathTracer): Task 6
    // Use Snell's Law to refract out_dir through the surface
    // Return the refracted direction. Set was_internal to false if
    // refraction does not occur due to total internal reflection,
    // and true otherwise.

    // When dot(out_dir,normal=(0,1,0)) is positive, then out_dir corresponds to a
    // ray exiting the surface into vaccum (ior = 1). However, note that
    // you should actually treat this case as _entering_ the surface, because
    // you want to compute the 'input' direction that would cause this output,
    // and to do so you can simply find the direction that out_dir would refract
    // _to_, as refraction is symmetric.
    float eta_t, eta_i;
    if (out_dir.y > 0)
        eta_t = index_of_refraction, eta_i = 1;
    else
        eta_t = 1, eta_i = index_of_refraction;
    float out_x = - out_dir.x * eta_i / eta_t, out_z = - out_dir.z * eta_i / eta_t;
    float out_y_sq = 1 - out_x * out_x - out_z * out_z;
    if (out_y_sq > 0) {
        was_internal = false;
        return Vec3(out_x, out_dir.y > 0 ? -sqrt(out_y_sq) : sqrt(out_y_sq), out_z);
    }
    else {
        was_internal = true;
        return out_dir;
    }
}

BSDF_Sample BSDF_Lambertian::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 5
    // Implement lambertian BSDF. Use of BSDF_Lambertian::sampler may be useful

    BSDF_Sample ret;
    ret.attenuation = out_dir.y > 0 ? albedo / PI_F : Spectrum(); // What is the ratio of reflected/incoming light?
    ret.direction = sampler.sample(ret.pdf);
    return ret;
}

Spectrum BSDF_Lambertian::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    return albedo * (1.0f / PI_F);
}

BSDF_Sample BSDF_Mirror::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement mirror BSDF

    BSDF_Sample ret;
    ret.attenuation = out_dir.y > 0 ? Spectrum(1.0) : Spectrum();
    ret.direction = reflect(out_dir);
    ret.pdf = 1.0f;
    return ret;
}

Spectrum BSDF_Mirror::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // Technically, we would return the proper reflectance
    // if in_dir was the perfectly reflected out_dir, but given
    // that we assume these are single exact directions in a
    // continuous space, just assume that we never hit them
    // _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Glass::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6

    // Implement glass BSDF.
    // (1) Compute Fresnel coefficient. Tip: use Schlick's approximation.
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?

    BSDF_Sample ret;
    bool was_internal;
    Vec3 r_dir = refract(out_dir, index_of_refraction, was_internal);
    ret.attenuation = Spectrum(1.0);
    ret.pdf = 1.0f;
    if (was_internal) {
        ret.direction = reflect(out_dir);
    } else {
        float eta_t, eta_i;
        if (out_dir.y > 0)
            eta_t = index_of_refraction, eta_i = 1;
        else
            eta_t = 1, eta_i = index_of_refraction;
        float s_i = std::abs(out_dir.y), s_t = std::abs(r_dir.y);
        float r_ll  = (eta_t * s_i - eta_i * s_t) / (eta_t * s_i + eta_i * s_t);
        float r__l_ = (eta_i * s_i - eta_t * s_t) / (eta_i * s_i + eta_t * s_t);
        float F_r = (r_ll * r_ll + r__l_ * r__l_) / 2;
        if (RNG::coin_flip(F_r))
            ret.direction = reflect(out_dir);
        else
            ret.direction = r_dir;
    }
    return ret;
}

Spectrum BSDF_Glass::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Diffuse::sample(Vec3 out_dir) const {
    BSDF_Sample ret;
    ret.direction = sampler.sample(ret.pdf);
    ret.emissive = radiance;
    ret.attenuation = {};
    return ret;
}

Spectrum BSDF_Diffuse::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // No incoming light is reflected; only emitted
    return {};
}

BSDF_Sample BSDF_Refract::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement pure refraction BSDF.

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?

    BSDF_Sample ret;
    bool was_internal;
    ret.direction = refract(out_dir, index_of_refraction, was_internal);
    ret.attenuation = was_internal? Spectrum() : Spectrum(1.0); // What is the ratio of reflected/incoming light?
    ret.pdf = 1.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)
    return ret;
}

Spectrum BSDF_Refract::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

} // namespace PT
