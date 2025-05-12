#ifndef HELLO_H
#define HELLO_H

#include "sim_random.h"
#include "image.h"
#include "params.h"
#include "scene.h"

typedef void render(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *random);

render render_image_phong;
render render_image_tracing;
render finite_differences;
render finite_differences_tracing;
render render_image_effects;

#endif // HELLO_H
