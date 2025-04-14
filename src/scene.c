#include <math.h>
#include "linmath.h"
#include "sdf.c"

typedef struct {
    float offset;
} SceneParams;

void params_from_float_pointer(const float *params, SceneParams *out) {
    out->offset = params[0];
}

const float *float_pointer_from_params(const SceneParams *out) {
    return &out->offset;
}

typedef struct {
    float distance;
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float shininess;
} SceneSample;

static void default_scene_sample(SceneSample *s) {
    s->distance = INFINITY;
    vec3_set(s->ambient, 0.0);
    vec3_set(s->diffuse, 0.0);
    vec3_set(s->specular, 0.0);
    s->shininess = 1.0;
}

/** writes the output in the first paramemter */
static void compose_scene_sample(SceneSample *destination, SceneSample *b) {
    if (destination->distance < b->distance) {
        return; // nothing to be done
    }
    destination->distance = b->distance;
    vec3_dup(destination->ambient, b->ambient);
    vec3_dup(destination->diffuse, b->diffuse);
    vec3_dup(destination->specular, b->specular);
    destination->shininess = b->shininess;
}

typedef void SDF_OBJECT(const vec3 pos, const SceneParams *params, SceneSample *sample);

void object_foreground_capsule(const vec3 pos, const SceneParams *params, SceneSample *sample) {
    vec3 offset = {-1.0, 1.0, 6.0};
    offset[0] += params->offset;
    vec3 pos1;
    vec3_add(pos1, pos, offset);
    sample->distance = sdfCylinder(pos1, 1.0f, 2.0f);
    vec3_set(sample->ambient, 0.1f);
    vec3 cDiffuse = {0.3f, 0.5f, 0.8f};
    vec3_dup(sample->diffuse, cDiffuse);
}


void object_foreground(const vec3 pos, const SceneParams *params, SceneSample *sample) {
    sample->distance = sdfVerticalCapsule(pos, 2.0, 1.0);
    vec3_set(sample->diffuse, 1.0f);
    vec3_set(sample->specular, 0.1f);
}

void scene_sample(const vec3 pos, const SceneParams *params, SceneSample *sample) {
    default_scene_sample(sample);

    SceneSample working;

    default_scene_sample(&working);
    object_foreground_capsule(pos, params, &working);
    compose_scene_sample(sample, &working);

    default_scene_sample(&working);
    object_foreground(pos, params, &working);
    compose_scene_sample(sample, &working);
}

//sdfCylinder(pos, 0.5, 1.0, param, result);
//vec3 normal = {1.0, 0.0, 0.0};
//sdfPlane(pos, normal, 0.3, 0.5, result);
// sdfVerticalCapsule(pos, 2.0, 1.0, param, result);
//sdfTriPrism(pos, 1.0, 5.0, param, result);
