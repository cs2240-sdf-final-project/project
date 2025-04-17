#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include "hello.h"

int main(void) {
    mkdir("debug-gradient", 0777); // do nothing if debug-gradient already exists

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);

    SceneParams *params = make_scene_params();

    render_image(&real, &gradient, params);

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

    free_image(&real);
    free_gradient_image(&gradient);
    free_scene_params(params);
}
