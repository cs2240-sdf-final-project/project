#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <assert.h>
#include <cctype>
#include <math.h>
#include "linmath.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include <functional>

#include "params.h"
#include "render.h"
#include "scene.h"
#include "util.h"
#include "sim_random.h"

#include "scene.cpp" // yes, enzyme needs to see these as part of the same translation unit

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// blur kernel effect radius
const int blur_radius = 5;

// we take steps at most this size in order to avoid missing
// sign changes in the directional derivatives
const float max_step_size = 1.0f;
// take this many steps to balance performance with exploring the entire scene
const int number_of_steps = 1'000;
// if our sdf gets smaller than this amount, we will consider it an intersection with the surface
const float contact_threshold = 1e-4f;
// width of the scene band
const float distance_threshold = 5e-2f;

const float lm_pi = 3.14159265358979323846f;

const float pathContinuationProb = 0.9f;
// finite differences half epsilon
const float finite_difference_epsilon = 1e-3f;

// how many samples per pixel in path tracing
const int numberOfSampling = 1;

enum BisectAffinity {
    BISECT_LEFT,
    BISECT_RIGHT,
    BISECT_STOP,
};

typedef BisectAffinity AFFINITY_FUNC(float evaluate_at, int iter_count, void *context);
typedef float MIDPOINT_FUNC(float a, float b, void *context);
inline float bisect(float t_min, float t_max,
    AFFINITY_FUNC *get_affinity, MIDPOINT_FUNC *get_midpoint, void *context) {
    int iter_count = 0;
    for (;;) {
        float t_mid = get_midpoint(t_min, t_max, context);
        BisectAffinity aff = get_affinity(t_mid, iter_count, context);
        if (aff == BISECT_STOP) {
            return t_mid;
        } else if (aff == BISECT_LEFT) {
            t_max = t_mid;
        } else if (aff == BISECT_RIGHT) {
            t_min = t_mid;
        }
        iter_count += 1;
    }
}

typedef struct {
    vec3 origin, direction;
    const SceneParams *params;
    const SceneContext *ctx;
    int iter_cap;
} CriticalPointContext;

inline float critical_point_bisect_midpoint(float a, float b, void *context) {
    CriticalPointContext *ctx = (CriticalPointContext *)context;
    float dir_a = directional_derivative(ctx->origin, ctx->direction, a, ctx->params, ctx->ctx);
    float dir_b = directional_derivative(ctx->origin, ctx->direction, b, ctx->params, ctx->ctx);
    // have points (a, dir_a) and (b, dir_b). want to find the zero crossing.
    float denom = dir_a - dir_b;
    if (fabsf(denom) < 1e-5) {
        return a; // arbitrarily choose one of the endpoints
    }
    float numer = dir_a * b - a * dir_b;
    return numer / denom;
}

inline BisectAffinity critical_point_bisect_affinity(float t, int iter_count, void *context) {
    CriticalPointContext *ctx = (CriticalPointContext *)context;
    if (iter_count >= ctx->iter_cap) {
        return BISECT_STOP;
    }
    // do a distance check on a fixed small iteration number to stop early
    const int early_stop_check_iter = 2;
    const float preliminary_distance_threshold = 1.0f;
    if (iter_count == early_stop_check_iter) {
        vec3 pos;
        ray_step(pos, ctx->origin, ctx->direction, t);
        if (scene_sample_sdf(pos, ctx->params, ctx->ctx) > preliminary_distance_threshold) {
            return BISECT_STOP;
        }
    }
    float dir_t = directional_derivative(ctx->origin, ctx->direction, t, ctx->params, ctx->ctx);
    bool is_approaching = dir_t < 0;
    return is_approaching ? BISECT_RIGHT : BISECT_LEFT;
}

typedef struct {
    bool found_critical_point;
    float t_if_found_critical_point;
} SearchResult;

void search_for_critical_point(
    const vec3 origin, const vec3 direction, const SceneParams *params, const SceneContext *ctx, float t_min, float t_max, SearchResult *ret
) {
    ret->found_critical_point = false;
    ret->t_if_found_critical_point = 0.0;

    // we consider directional derivatives this large to be zero
    const float directional_derivative_threshold = 1e-1f;
    // we will bisect this many iterations in order to find the true location of the minimum directional derivative
    const int extensive_depth = 12;
    CriticalPointContext context;
    vec3_dup(context.origin, origin);
    vec3_dup(context.direction, direction);
    context.params = params;
    context.iter_cap = extensive_depth;
    context.ctx = ctx;
    float best_t = bisect(t_min, t_max, critical_point_bisect_affinity, critical_point_bisect_midpoint, &context);
    vec3 best_pos;
    ray_step(best_pos, origin, direction, best_t);
    if (scene_sample_sdf(best_pos, params, ctx) > distance_threshold) {
        return;
    }
    float best_dir_t = directional_derivative(origin, direction, best_t, params, ctx);
    if (fabsf(best_dir_t) > directional_derivative_threshold) {
        return;
    }
    ret->found_critical_point = true;
    ret->t_if_found_critical_point = best_t;
}

float sdf_theta_wrapper(const vec3 pos, const float *params, const SceneContext *ctx) {
    SceneParams *scene_params = params_from_float_pointer(params);
    return scene_sample_sdf(pos, scene_params, ctx);
}

extern void __enzyme_autodiff_theta(void *, int,const float *, int, const float *, float *, int, const SceneContext *ctx);

void diff_sdf(const vec3 pos, SceneParams *paramsOut, const SceneParams *paramsIn, const SceneContext *ctx) {
    const float *raw_params = float_pointer_from_params(paramsIn);
    float* draw_params = float_pointer_from_params(paramsOut);
    __enzyme_autodiff_theta(
        (void*)sdf_theta_wrapper,
        enzyme_const, pos,
        enzyme_dup, raw_params, draw_params,
        enzyme_const, ctx
    );
}

inline void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SdfResult *sample, const SceneParams *params) {
    (void)params;
    float lightColors[3][3] = {
        {.2f, .2f, .2f},
        {.2f, .2f, .2f},
        {.2f, .2f, .2f},
    };
    // vec3_add(lightColors[0], lightColors[0], params->color);
    // world[0] -= params->object_3_pos[0];
    float lightDirections[3][3] = {
        {0.f, -1.f, 0.f},
        {-3.f, 2.f, 0.f},
        {0.f, -1.f, 3.f},
    };
    float kd = 0.5;
    float ks = 1.0;
    vec3_dup(radiance, sample->ambient);
    for (int l = 0; l < 3; l++) {
        vec3 lightColor;
        vec3_dup(lightColor, lightColors[l]);
        // lightColor[0] += params->object_1_color[0];
        // lightColor[1] += params->object_1_color[1];
        // lightColor[2] += params->object_1_color[2];
        vec3 light_dir;
        vec3_norm(light_dir, lightDirections[l]);
        // light_dir[0] += params->object_1_direction[0];
        // light_dir[1] += params->object_1_direction[1];
        // light_dir[2] += params->object_1_direction[2];
        float facing = fmaxf(0.0, vec3_mul_inner(normal, light_dir));


        vec3 bounce;
        for (int i = 0; i < 3; i++) {
            bounce[i] = light_dir[i] - normal[i] * facing * 2.f;
        }
        float specular = powf(vec3_mul_inner(bounce, looking), sample->shininess);
        for (int i = 0; i < 3; i++) {
            radiance[i] += kd * facing * sample->diffuse[i] * lightColor[i];
            radiance[i] += ks * specular * sample->specular[i] * lightColor[i];
        }
    }
}

typedef struct {
    bool found_intersection;
    float intersection_t;
} IntersectionResult;

inline IntersectionResult trace_ray_get_intersection(
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    float t = 0.0;
    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params, ctx);
        if (distance < contact_threshold) {
            ret.found_intersection = true;
            ret.intersection_t = t;
            break;
        }
        t += distance;
    }
    return ret;
}

IntersectionResult trace_ray_get_critical_point(
    SearchResult *critical_point,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    float t = 0.0;
    float previous_t = 0.0;
    critical_point->found_critical_point = false;
    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        float dir_previous_t = directional_derivative(origin, direction, previous_t, params, ctx);
        float dir_t = directional_derivative(origin, direction, t, params, ctx);
        // look for a sign change in the directional derivative
        if (!critical_point->found_critical_point && ((dir_previous_t < 0) && (dir_t > 0))) {
            // let's try to find critical_point between t and previous_t;
            search_for_critical_point(origin, direction, params, ctx, previous_t, t, critical_point);
        }
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params, ctx);
        if (distance < contact_threshold) {
            ret.found_intersection = true;
            ret.intersection_t = t;
            break;
        }
        float step_size = fminf(max_step_size, distance);
        previous_t = t;
        t += step_size;
    }
    return ret;
}

float normal_pdf(float x, float mu, float sigma) {
    if (sigma <= 0.0f) { return 0.0f; } // Avoid division by zero / invalid input
    float diff = x - mu;
    float exponent = -(diff * diff) / (2.0f * sigma * sigma);
    float scaling_factor = 1.0f / (sigma * sqrtf(2.0f * lm_pi));
    return scaling_factor * expf(exponent);
}

// --------------------------------------------------------------------

void uniformSampleBottomFaceOriginalSDF(int num_points, vec3 *output_points, RandomState* random) {
    for (int i = 0; i < num_points; ++i) {
        float randomX = random_next_float(random);
        float randomZ = random_next_float(random);
        float sampledX = -0.5f + randomX * 1.0f;
        float sampledZ = -1.f * 0.5f + randomZ * 1.f;
        vec3 result;
        result[0] = sampledX;
        result[1] = -1;
        result[2] = sampledZ;
        vec3_dup(output_points[i],result);

    }
}

void computeDirectLighting(vec3 hitPoint, vec3 hitNormal, SdfResult sample, const SceneParams *params, const SceneContext *ctx, vec3 directLighting, RandomState* random){
    const int num_samples = 1;
    vec3 sampled_points[num_samples];

    vec3 center;
    vec3_set(center, 0.0f, 1.f, 0.0f); // Center of distribution on the plane y=1.0

    uniformSampleBottomFaceOriginalSDF(num_samples,sampled_points,random);

    vec3 contribution;
    vec3_set(contribution, 0.0f);

    for(int i = 0; i < num_samples; i++){
        vec3 currentSamplePoint;

        vec3_dup(currentSamplePoint,sampled_points[i]);

        vec3 light_direction;

        vec3_sub(light_direction, currentSamplePoint, hitPoint);

        float distance2 = vec3_mul_inner(light_direction, light_direction);

        vec3_norm(light_direction, light_direction); // point to light

        IntersectionResult shadow_hit_point = trace_ray_get_intersection(hitPoint, light_direction, params, ctx);
        if (!shadow_hit_point.found_intersection) {
            // trace_ray_get_intersection(hitPoint, light_direction, params, ctx, true);
        }
        if (shadow_hit_point.found_intersection) {
            vec3 shadow_hit_intersection_point;
            ray_step(shadow_hit_intersection_point, hitPoint, light_direction, shadow_hit_point.intersection_t);

            vec3 shawdow_normal;
            get_normal_from(shawdow_normal,shadow_hit_intersection_point,params, ctx);

            vec3 diff_hit;
            vec3_sub(diff_hit,shadow_hit_intersection_point,hitPoint);
            float hitDistance = vec3_len(diff_hit);

            vec3 diff_light;
            vec3_sub(diff_light,currentSamplePoint,hitPoint);
            float distToLight = vec3_len(diff_light);

            if (fabsf(hitDistance - distToLight) < 1e-4f) { // not occluded...

                float flip = vec3_mul_inner(light_direction,hitNormal);

                if (flip < 0) {
                    vec3_scale(hitNormal,hitNormal,-1.0f);
                }

                float costheta = fmaxf(0.f, vec3_mul_inner(light_direction,hitNormal));

                vec3 emission;
                vec3_set(emission, 10.f);
                vec3 lightNormal;
                vec3_set(lightNormal, 0.0f, 1.f, 0.0f);

                float cosThetaPrime = fmax(0.0f, -vec3_mul_inner(light_direction, lightNormal));

                vec3 bsdf;
                vec3_set(bsdf,1.0f);

                vec3_scale(bsdf, sample.diffuse, 1.0f / lm_pi);

                float lightArea = 1.0f * 0.8f;
                float pdf = 1.0f / lightArea;

                // theta
                float last_term = costheta * cosThetaPrime / ((distance2 + 1e-6f) * pdf);

                vec3_scale(bsdf, bsdf, last_term);
                vec3_cwiseProduct(bsdf, bsdf, emission);

                vec3_add(contribution, contribution, bsdf);
            }
        }
        else {
            // vec3 blue = {0.f, 0.f, 100.f};
            // vec3_add(contribution, contribution, blue);
        }
    }
    vec3_scale(contribution, contribution, 1.0f/num_samples);
    vec3_dup(directLighting, contribution);
}

void gammaCorrect(vec3 gamma_correct_color, const vec3 color) {
    float gamma = 2.2f;
    gamma_correct_color[0] = powf(fmaxf(0.0f, color[0]), 1.0f / gamma);
    gamma_correct_color[1] = powf(fmaxf(0.0f, color[1]), 1.0f / gamma);
    gamma_correct_color[2] = powf(fmaxf(0.0f, color[2]), 1.0f / gamma);
}

inline void get_radiance_phong(
    vec3 radiance,
    const IntersectionResult *intersection,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    if (!intersection->found_intersection) {
        vec3_set(radiance, 0.f);
        return;
    }
    vec3 current_position;
    ray_step(current_position, origin, direction, intersection->intersection_t);
    SdfResult sample;
    scene_sample(current_position, params, ctx, &sample);
    vec3 normal;
    get_normal_from(normal, current_position, params, ctx);
    phongLight(radiance, direction, normal, &sample, params);
}

float render_get_radiance_wrapper(
    const IntersectionResult *intersection,
    const vec3 origin,
    const vec3 direction,
    const float *raw_params,
    const SceneContext *ctx,
    int ch
) {
    const SceneParams *params = params_from_float_pointer(raw_params);
    vec3 naive_intersection;
    ray_step(naive_intersection, origin, direction, intersection->intersection_t);
    float sdf_at_intersection = scene_sample_sdf(naive_intersection, params, ctx);
    float dsdfdt = directional_derivative(origin, direction, intersection->intersection_t, params, ctx);
    float first_order_update = - sdf_at_intersection / dsdfdt;
    IntersectionResult sensitivity_adjusted;
    sensitivity_adjusted.found_intersection = intersection->found_intersection;
    sensitivity_adjusted.intersection_t = intersection->intersection_t + first_order_update;
    vec3 radiance;
    get_radiance_phong(radiance, &sensitivity_adjusted, origin, direction, params, ctx);
    return radiance[ch];
}

void render_pixel_phong(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    IntersectionResult intersection = trace_ray_get_intersection(origin, direction, params, ctx);
    get_radiance_phong(real, &intersection, origin, direction, params, ctx);
}

extern void __enzyme_autodiff_radiance(
    void *,
    int, const IntersectionResult *,
    int, const float *,
    int, const float *,
    int, const float *, float *,
    int, const SceneContext *ctx,
    int, int
);

typedef struct {
    vec3 pos;
    vec3 dir;
    vec3 contribution;
} Segment;

void step_segment(Segment &segment, float t) {
    ray_step(segment.pos, segment.pos, segment.dir, t);
}

std::vector<Segment> getSecondaryPath(
    const vec3 origin,
    const vec3 direction,
    RandomState *random,
    const SceneParams *params,
    const SceneContext *ctx
) {
    std::vector<Segment> path;

    vec3 current_position;
    vec3_dup(current_position, origin);
    vec3 current_direction;
    vec3_dup(current_direction,direction);
    vec3 previous_direction;
    vec3_set(previous_direction,0.0f);

    int iter = 0;

    // While the ray has not terminated...
    while(true){
        iter++;


        // Make a new segment from origin in direction of hit point
        Segment newSegment;
        vec3_dup(newSegment.pos, current_position);
        vec3_dup(newSegment.dir, current_direction);

        vec3 directLighting; // it re-declares directLighting every iteration/every bounce, is this ok?

        if(iter == 0){ // this will never happen...
            vec3_set(directLighting, 0.0f);
        }
        else{
            SdfResult hit_sample;
            scene_sample(current_position, params, ctx, &hit_sample);

            vec3 hit_normal;
            get_normal_from(hit_normal, current_position, params, ctx);
            computeDirectLighting(current_position,hit_normal,hit_sample,params,ctx, directLighting,random);

        }

        vec3_dup(newSegment.contribution, directLighting);

        path.push_back(newSegment);

        if(random_next_float(random)> pathContinuationProb){
            break;
        }

        // Find out if the ray hit
        IntersectionResult hitPoint = trace_ray_get_intersection(current_position,current_direction,params, ctx);
        if(!hitPoint.found_intersection){
            return path;
        }

        //find intersection point

        vec3 hit_intersection_point;
        ray_step(hit_intersection_point, current_position, current_direction, hitPoint.intersection_t);

        SdfResult sample;
        scene_sample(hit_intersection_point, params, ctx, &sample);

        vec3 normal;
        get_normal_from(normal, hit_intersection_point, params, ctx);


        vec3 newDir;
        if(sample.isReflected){
            vec3 refl_dir;
            vec3_reflect(refl_dir,current_direction,normal);
            vec3_dup(newDir,refl_dir);
            //sample_hemisphere(newDir, normal, random);
        }else{

            //uniformSampleHemisphere(newDir, normal);
            sample_hemisphere(newDir, normal, random);

        }
        //sample_hemisphere(newDir, normal, random);
        //uniformSampleHemisphere(newDir, normal);

        const float ray_epsilon = 1e-4f;

        vec3 offset_position;
        vec3_scale(offset_position, normal, ray_epsilon);
        vec3_add(hit_intersection_point,hit_intersection_point,offset_position);

        //update new ray
        vec3_dup(current_position,hit_intersection_point);
        vec3_dup(previous_direction,current_direction);
        vec3_dup(current_direction,newDir);
    }

    // std::cout << std::endl;

    return path;
}

void get_radiance_tracing(
    vec3 radiance,
    const std::vector<Segment> &path,
    const SceneParams *params,
    const SceneContext *ctx
) {
    vec3 spectralFilter;
    vec3_set(spectralFilter, 1.f);
    vec3 intensity;
    vec3_set(intensity, 0.f);
    for(size_t i = 0; i < path.size(); i++){
        if(i != path.size() - 1){
            vec3 hitPosition;
            vec3_dup(hitPosition, path[i+1].pos);
            vec3 hitDirection;
            vec3_dup(hitDirection,path[i].dir);
            vec3 wi;
            vec3_dup(wi,path[i+1].dir);

            SdfResult sample;
            scene_sample(hitPosition, params, ctx, &sample);
            vec3 normal;
            get_normal_from(normal, hitPosition, params, ctx);

            vec3 emissive;
            vec3_dup(emissive, sample.emissive);

            vec3 emissive_part;
            vec3_cwiseProduct(emissive_part, emissive, spectralFilter);

            vec3_add(intensity, intensity, emissive_part);

            if(sample.isReflected){
                vec3_scale(spectralFilter,spectralFilter, 1.0f / pathContinuationProb);
                vec3 brdf;
                vec3_set(brdf, 1.0f);
                vec3_cwiseProduct(spectralFilter, spectralFilter, brdf);
            }else{
                vec3_scale(spectralFilter,spectralFilter, 1.0f / pathContinuationProb);

                float pdf = 2.f * lm_pi;
                vec3_scale(spectralFilter,spectralFilter, pdf);

                float cosine_term = fmaxf(0.f, vec3_mul_inner(normal, wi));
                    //float cosine_term = vec3_mul_inner(normal, wi);
                vec3_scale(spectralFilter, spectralFilter, cosine_term);

                vec3 brdf;
                vec3_scale(brdf, sample.diffuse, 1.0f / lm_pi);
                vec3_cwiseProduct(spectralFilter, spectralFilter, brdf);
            }
        }
    }
    vec3_dup(radiance,intensity);
}

float render_pixel_tracing_wrapper(
    const std::vector<Segment> &path,
    const float *raw_params,
    const SceneContext *ctx,
    int ch
) {
    const SceneParams *params = params_from_float_pointer(raw_params);
    vec3 radiance;
    get_radiance_tracing(radiance, path, params, ctx);
    return radiance[ch];
}

void render_gradient_phong(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    const SceneContext *ctx,
    RandomState *rng
) {
    (void)rng;
    SearchResult critical_point;
    IntersectionResult intersection = trace_ray_get_critical_point(&critical_point, origin, direction, params, ctx);
    const float *raw_params = float_pointer_from_params(params);
    fill_ppc(params_per_channel, 0.f);
    for (int ch = 0; ch < 3; ch++) {
        float *raw_fill = ppc_get_channel(params_per_channel, ch);
        extern void __enzyme_autodiff_radiance(
            void *,
            int, const IntersectionResult *,
            int, const float *,
            int, const float *,
            int, const float *, float *,
            int, const SceneContext *ctx,
            int, int
        );
        __enzyme_autodiff_radiance(
            (void*)render_get_radiance_wrapper,
            enzyme_const, &intersection,
            enzyme_const, origin,
            enzyme_const, direction,
            enzyme_dupnoneed, raw_params, raw_fill,
            enzyme_const, ctx,
            enzyme_const, ch
        );
    }
    get_radiance_phong(real, &intersection, origin, direction, params, ctx);
    if (critical_point.found_critical_point) {
        SceneParams *dummy_params = uninit_scene_params();
        scene_params_fill(dummy_params, 0.f);
        vec3 y_star;
        ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
        SdfResult sample;
        scene_sample(y_star, params, ctx, &sample);
        vec3 normal;
        get_normal_from(normal, y_star, params, ctx);
        vec3 y_star_radiance;
        phongLight(y_star_radiance, direction, normal, &sample, params);
        diff_sdf(y_star, dummy_params, params, ctx);
        vec3 deltaL;
        vec3_sub(deltaL, y_star_radiance, real);
        float extra_kick = 50.f;
        vec3_scale(deltaL, deltaL, -1.f * extra_kick / distance_threshold);
        outer_product_add_assign(params_per_channel, dummy_params, deltaL);
        free_scene_params(dummy_params);
    }
    params_nan_to_num(params_per_channel, 0.0);
}

void render_pixel_tracing(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    const SceneContext *ctx,
    RandomState *rng
) {
    SearchResult critical_point;
    IntersectionResult intersection = trace_ray_get_critical_point(&critical_point, origin, direction, params, ctx);
    (void)intersection;
    std::vector<Segment> path = getSecondaryPath(origin, direction, rng, params, ctx);
    get_radiance_tracing(real, path, params, ctx);
    const float *raw_params = float_pointer_from_params(params);
    fill_ppc(params_per_channel, 0.0);
    for (int ch = 0; ch < 3; ch++) {
        float *raw_fill = ppc_get_channel(params_per_channel, ch);
        extern void __enzyme_autodiff_radiance2(
            void *,
            int, const std::vector<Segment> &,
            int, const float *, float *,
            int, const SceneContext *ctx,
            int, int
        );
        __enzyme_autodiff_radiance2(
            (void*)render_pixel_tracing_wrapper,
            enzyme_const, path,
            enzyme_dupnoneed, raw_params, raw_fill,
            enzyme_const, ctx,
            enzyme_const, ch
        );
    }
    if (critical_point.found_critical_point && path.size() > 0) {
        SceneParams *dummy_params = uninit_scene_params();
        scene_params_fill(dummy_params, 0.f);
        vec3 y_star;
        ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
        vec3 y_star_radiance;
        std::vector<Segment> other_path(path);
        step_segment(other_path.at(0), critical_point.t_if_found_critical_point);
        get_radiance_tracing(y_star_radiance, other_path, params, ctx);
        diff_sdf(y_star, dummy_params, params, ctx);
        vec3 deltaL;
        vec3_sub(deltaL, y_star_radiance, real);
        vec3_scale(deltaL, deltaL, -1.f / distance_threshold);
        outer_product_add_assign(params_per_channel, dummy_params, deltaL);
        free_scene_params(dummy_params);
    }
    params_nan_to_num(params_per_channel, 0.0);
}

void render_pixel_tracing_wrapper_super_sample(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    const SceneContext *ctx,
    RandomState *rng
) {
    SceneParamsPerChannel* working = make_ppc();

    vec3_set(real, 0.0);
    fill_ppc(params_per_channel, 0.0);
    for(int i = 0; i < numberOfSampling; i++){
        vec3 temp_real;
        render_pixel_tracing(temp_real, origin,direction, working, params, ctx,rng);
        vec3_add(real,real, temp_real);
        params_per_channel_add_assign(params_per_channel, working);
    }
    vec3_scale(real,real, 1.f/numberOfSampling);
    free_ppc(working);
}

void radiance_tracing_super_sample(vec3 real,
        const vec3 origin,
        const vec3 direction,
        const SceneParams *params,
        const SceneContext *ctx,
        RandomState *rng){

    vec3_set(real, 0.0);

    for(int i = 0; i < numberOfSampling; i++){
        vec3 temp_real;
        std::vector<Segment> path_left = getSecondaryPath(origin, direction, rng, params, ctx);
        get_radiance_tracing(temp_real,path_left,params,ctx);
        vec3_add(real,real, temp_real);
    }
    vec3_scale(real,real, 1.f/numberOfSampling);
}

struct Renderer {
    float world_from_camera[4][4];
    vec3 camera_position;
    long image_width, image_height;
};

Renderer *make_renderer(long image_width, long image_height) {
    float aspect = (float)image_width / (float)image_height ;
    float near_clip = 0.1f;
    float far_clip = 100.0f;
    float y_fov = 0.785f;
    vec3 camera_position = {0.f, 0.f, 3.6f};
    vec3 center = {0};
    vec3 up = {0.0, 1.0, 0.0};
    mat4x4 look_at;
    mat4x4_look_at(look_at, camera_position, center, up);
    mat4x4 perspective;
    mat4x4_perspective(perspective, y_fov, aspect, near_clip, far_clip);
    mat4x4 projection_view;
    mat4x4_mul(projection_view, perspective, look_at);
    mat4x4 world_from_camera;
    mat4x4_invert(world_from_camera, projection_view);
    Renderer *ret = new Renderer;
    mat4x4_dup(ret->world_from_camera, world_from_camera);
    vec3_dup(ret->camera_position, camera_position);
    ret->image_width = image_width;
    ret->image_height = image_height;
    return ret;
}

void free_renderer(Renderer *renderer) {
    delete renderer;
}

struct PixelRenderer {
    float camera_position[3];
    float direction[3];
};

PixelRenderer make_pixel_renderer(const Renderer *renderer, long ir, long ic) {
    float r = (float)ir;
    float c = (float)ic;
    float device_x = lerp(c, 0.0, (float)renderer->image_width, -1.0, 1.0);
    float device_y = lerp(r, 0.0, (float)renderer->image_height, -1.0, 1.0);
    vec4 unprojected = {device_x, device_y, 1.0, 1.0};
    vec4 homo;
    mat4x4_mul_vec4(homo, renderer->world_from_camera, unprojected);
    vec3 far;
    dehomogenize(far, homo);
    vec3 direction;
    vec3_sub(direction, far, renderer->camera_position);
    vec3_norm(direction, direction);
    PixelRenderer ret;
    vec3_dup(ret.camera_position, renderer->camera_position);
    vec3_dup(ret.direction, direction);
    return ret;
}

void render_image_phong(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;

    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);

    Renderer *renderer = make_renderer(image_width, image_height);
    #pragma omp parallel for
    for (long ir = 0; ir < image_height; ir++) {
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            vec3 out_real;
            render_gradient_phong(out_real, pr.camera_position, pr.direction, ppc, params, ctx, rng);
            image_set(real, ir, ic, out_real);
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void render_image_tracing(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;
    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    Renderer *renderer = make_renderer(image_width, image_height);
    #pragma omp parallel for
    for (long ir = 0; ir < image_height; ir++) {
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            vec3 out_real;
            render_pixel_tracing(out_real, pr.camera_position, pr.direction, ppc, params, ctx, rng);
            image_set(real, ir, ic, out_real);
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void free_pixel_renderer(PixelRenderer *renderer) {
    delete renderer;
}

void finite_differences(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;
    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    Renderer *renderer = make_renderer(image_width, image_height);

    #pragma omp parallel for
    for (long ir = 0; ir < image_height; ir++) {
        SceneParams *working = uninit_scene_params();
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            for (long p = 0; p < number_of_scene_params; p++) {
                float param_value = scene_parameter_get(params, p);
                scene_params_copy(working, params);
                scene_params_set(working, p, param_value - finite_difference_epsilon);
                vec3 real_left;
                render_pixel_phong(real_left, pr.camera_position, pr.direction, working, ctx);
                vec3 real_right;
                scene_params_copy(working, params);
                scene_params_set(working, p, param_value + finite_difference_epsilon);
                render_pixel_phong(real_right, pr.camera_position, pr.direction, working, ctx);

                vec3 gradient;
                vec3_sub(gradient, real_right, real_left);
                vec3_scale(gradient, gradient, 1.f / (2.f * finite_difference_epsilon));

                set_ppc(ppc, p, gradient);
            }
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_scene_params(working);
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void finite_differences_tracing(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng){
    (void)rng;
    const float half_epsilon = 1e-1f;

    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    Renderer *renderer = make_renderer(image_width, image_height);
    #pragma omp parallel for
    for (long ir = 0; ir < image_height; ir++) {
        printf("current row is: %ld",ir);
        printf(" \n");
        fflush(stdout);
        SceneParams *working = uninit_scene_params();
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            for (long p = 0; p < number_of_scene_params; p++) {
                float param_value = scene_parameter_get(params, p);
                scene_params_copy(working, params);
                scene_params_set(working, p, param_value - half_epsilon);

                vec3 real_left;
                radiance_tracing_super_sample(real_left,pr.camera_position,pr.direction,working,ctx,rng);

                scene_params_copy(working, params);
                scene_params_set(working, p, param_value + half_epsilon);

                vec3 real_right;
                radiance_tracing_super_sample(real_right,pr.camera_position,pr.direction,working,ctx,rng);

                vec3 gradient;
                vec3_sub(gradient, real_right, real_left);
                vec3_scale(gradient, gradient, 1.f / (half_epsilon));
                set_ppc(ppc, p, gradient);
            }
            std::vector<Segment> path_debug = getSecondaryPath(pr.camera_position, pr.direction, rng, params, ctx);
            vec3 real_debug;
            get_radiance_tracing(real_debug,path_debug, params, ctx);
            image_set(real, ir, ic, real_debug);
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_ppc(ppc);
        free_scene_params(working);
    }
    free_renderer(renderer);
}

void chromaticAberrationFilter(long rShift, long gShift, long bShift, Image* outImage, const Image* real){
    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    for (long ir = 0; ir < image_height; ir++) {
        for (long ic = 0; ic < image_width; ic++) {
            vec3 pixel;
            image_get(pixel,real,ir,ic);

            long rShiftedCol = lclamp(ic + rShift, 0, image_width - 1);
            long gShiftedCol = lclamp(ic + gShift, 0, image_width - 1);
            long bShiftedCol = lclamp(ic + bShift, 0, image_width - 1);

            vec3 rChannelPixel;
            vec3 gChannelPixel;
            vec3 bChannelPixel;
            image_get(rChannelPixel, real, ir, rShiftedCol);
            image_get(gChannelPixel, real, ir, gShiftedCol);
            image_get(bChannelPixel, real, ir, bShiftedCol);

            vec3 outputPixel;
            outputPixel[0] = rChannelPixel[0];
            outputPixel[1] = gChannelPixel[1];
            outputPixel[2] = bChannelPixel[2];

            image_set(outImage, ir, ic, outputPixel);
        }
    }
}

float* generateBlurkernel(float sigma){
    float* kernel = nullptr;
    int kernel_size = 2 * blur_radius + 1;
    assert(kernel_size > 0);
    kernel = new float[static_cast<size_t>(kernel_size)];
    float sum = 0.0f;
    // Calculate the kernel values
    for (int i = 0; i < kernel_size; ++i) {
        float x = static_cast<float>(i - 2);  // Adjust x to the correct range
        float exponent = -(x * x) / (2 * sigma * sigma);
        kernel[i] = (1.0f / (sqrt(2 * lm_pi) * sigma)) * exp(exponent);
        sum += kernel[i];  // Accumulate the sum for normalization
    }
    for (int i = 0; i < kernel_size; ++i) {
        kernel[i] /= sum;
    }
    return kernel;
}

void blurFilter(Image* outImage, const Image* real, float sigma) {
    float* kernel = generateBlurkernel(sigma);
    long kernel_size = 2 * blur_radius + 1;
    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    for (long ir = 0; ir < image_height; ir++) {
        for (long ic = 0; ic < image_width; ic++) {
            float redAcc = 0.0f, greenAcc = 0.0f, blueAcc = 0.0f;
            for (long k = 0; k < kernel_size; ++k) {
                long pixelCol = lclamp(ic + k - blur_radius, 0L, static_cast<long>(image_width - 1));
                vec3 pixel;
                image_get(pixel,real,ir,pixelCol);
                float blurWeight = kernel[k];
                redAcc += pixel[0] * blurWeight;
                greenAcc += pixel[1]* blurWeight;
                blueAcc += pixel[2] * blurWeight;
            }
            vec3 out;
            out[0] = fmaxf(0.0f, fminf(1.f, redAcc));
            out[1] = fmaxf(0.0f, fminf(1.f, greenAcc));
            out[2] = fmaxf(0.0f, fminf(1.f, blueAcc));
            image_set(outImage, ir, ic, out);
        }
    }
    delete[] kernel;
}

void render_image_effects(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)gradient;
    long image_width, image_height;
    image_fetch_dims(real, &image_width, &image_height);
    Image *scratch = make_image(image_width, image_height);
    render_image_phong(real, gradient, params, ctx, rng);
    // chromaticAberrationFilter(5, 2, 10, &scratch, real);
    float sigma = 1.5f;
    blurFilter(scratch, real, sigma);
    image_copy(scratch, real);
    free_image(scratch);
}
