#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include "hello.h"

int main(void) {
    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);

    SceneParams *params = make_scene_params();
    params->offset = 0.1f;

    render_image(&real, &gradient, params);

    FILE *freal = fopen("real.bpm", "w");
    image_write_bpm(&real, freal);

    Image grad_slice = make_image(image_width, image_height);
    for (long p = 0; p < number_of_scene_params; p++) {
        std::ostringstream filename;
        filename << "gradient_" << std::setw(3) << std::setfill('0') << p << ".bpm";
        std::string filename_string = filename.str();
        FILE *fgradient = fopen(filename_string.c_str(), "w");
        gradient_image_slice(&grad_slice, &gradient, p);
        image_write_bpm(&grad_slice, fgradient);
    }

    free_image(&real);
    free_gradient_image(&gradient);
    free_scene_params(params);
}
