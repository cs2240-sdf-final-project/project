#include <fstream>
#include <iomanip>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <functional>
#include <sys/stat.h>
#include <filesystem>

#include <stdio.h>
#include <stdlib.h>
#include "hello.h"
#include "sim_random.h"

int main(void) {
    long image_width = 500;
    long image_height = 500;
    RandomState *rng = make_random();
    SceneContext *ctx = make_scene_context();
    SceneParams *params = uninit_scene_params();
    scene_params_init(params, ctx);
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);
    Image grad_slice = make_image(image_width, image_height);

    auto body = [&](std::filesystem::path dir, std::function<void(Image *, GradientImage *, const SceneParams *, const SceneContext *, RandomState *)> call) {

        call(&real, &gradient, params, ctx, rng);
        std::filesystem::create_directory(dir);
        std::filesystem::path file_path = dir / "real.ppm";
        FILE *freal = fopen(file_path.c_str(), "w");
        assert(freal);
        image_write_ppm(&real, freal);

        for (long p = 0; p < number_of_scene_params; p++) {
            std::ostringstream filename;
            filename << std::setw(3) << std::setfill('0') << p << ".ppm";
            std::filesystem::path grad_path = dir / filename.str();
            FILE *fgradient = fopen(grad_path.c_str(), "w");
            gradient_image_slice(&grad_slice, &gradient, p);
            image_write_ppm(&grad_slice, fgradient);
        }
    };

    body("debug-tracing", render_image_tracing);
    body("debug-phong", render_image_phong);
    body("debug-fd", finite_differences);
    body("debug-effects", render_image_effects);

    free_image(&grad_slice);
    free_image(&real);
    free_gradient_image(&gradient);
    free_random(rng);
    free_scene_params(params);
    free_scene_context(ctx);
}
