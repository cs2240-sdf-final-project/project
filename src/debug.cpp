#include <fstream>
#include <iomanip>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include "hello.h"
#include "sim_random.h"

static void debug_fd(const SceneParams *params, const SceneContext *ctx, RandomState *rng)  {
    mkdir("debug-fd", 0777); // do nothing if debug-fd already exists

    long image_width = 500;
    long image_height = 500;

    GradientImage gradient = make_gradient_image(image_width, image_height);
    finite_differences(&gradient, image_width, image_height, params, ctx, rng);

    Image grad_slice = make_image(image_width, image_height);
    for (long p = 0; p < number_of_scene_params; p++) {
        std::ostringstream filename;
        filename << "debug-fd/fd_" << std::setw(3) << std::setfill('0') << p << ".ppm";
        std::string filename_string = filename.str();
        FILE *fgradient = fopen(filename_string.c_str(), "w");
        gradient_image_slice(&grad_slice, &gradient, p);
        image_write_ppm(&grad_slice, fgradient);
    }

    free_image(&grad_slice);
    free_gradient_image(&gradient);
}

static void debug_gradient(const SceneParams *params, const SceneContext *ctx, RandomState *rng)  {
    mkdir("debug-gradient", 0777); // do nothing if debug-gradient already exists

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);

    GradientImage gradient = make_gradient_image(image_width, image_height);

    render_image(&real, &gradient, params, ctx, rng);

    FILE *freal = fopen("real.ppm", "w");
    image_write_ppm(&real, freal);

    Image grad_slice = make_image(image_width, image_height);
    for (long p = 0; p < number_of_scene_params; p++) {
        std::ostringstream filename;
        filename << "debug-gradient/gradient_" << std::setw(3) << std::setfill('0') << p << ".ppm";
        std::string filename_string = filename.str();
        FILE *fgradient = fopen(filename_string.c_str(), "w");
        gradient_image_slice(&grad_slice, &gradient, p);
        image_write_ppm(&grad_slice, fgradient);
    }

    free_image(&grad_slice);
    free_image(&real);
    free_gradient_image(&gradient);
}

int main(void) {
    RandomState *rng = make_random();

    SceneContext *ctx = make_scene_context();
    SceneParams *params = uninit_scene_params();
    scene_params_init(params, ctx);
    debug_gradient(params, ctx, rng);
    debug_fd(params, ctx, rng);

    free_random(rng);
    free_scene_params(params);
    free_scene_context(ctx);
}
