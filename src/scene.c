#include <math.h>
#include "scene.h"
#include "linmath.h"

typedef struct {
    float distance;
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float shininess;
} SceneSample;

static void fill_vec3(float *items, float fill_with) {
    for (int i = 0; i < 3; i++) {
        items[i] = fill_with;
    }
}

static void copy_vec3(float *destination, float *source) {
    for (int i = 0; i < 3; i++) {
        destination[i] = source[i];
    }
}

void default_scene_sample(SceneSample *s) {
    s->distance = INFINITY;
    fill_vec3(s->ambient, 0.0);
    fill_vec3(s->diffuse, 0.0);
    fill_vec3(s->specular, 0.0);
    s->shininess = 1.0;
}

/** writes the output in the first paramemter */
void compose_scene_sample(SceneSample *destination, SceneSample *b) {
    if (destination->distance > b->distance) {
        destination->distance = b->distance;
        copy_vec3(destination->ambient, b->ambient);
        copy_vec3(destination->diffuse, b->diffuse);
        copy_vec3(destination->specular, b->specular);
        destination->shininess = b->shininess;
    }
}
