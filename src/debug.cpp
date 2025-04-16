#include <stdio.h>
#include <stdlib.h>
#include "hello.h"

int main(void) {
    FILE *freal = fopen("real.bpm", "w");
    FILE *fgradient = fopen("gradient.bpm", "w");

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);


    SceneParams *params = make_scene_params();
    params->offset = 0.1f;

    render_image(&real, &gradient, params);
    image_write_bpm(&real, freal);
    // image_write_bpm(&gradient, fgradient);
    free_image(&real);
    // free_image(&gradient);
    free_scene_params(params);
}
