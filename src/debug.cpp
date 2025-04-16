#include <stdio.h>
#include <stdlib.h>
#include "hello.h"

int main(int argc, char *argv[]) {
    FILE *freal = fopen("real.bpm", "w");
    FILE *fgradient = fopen("gradient.bpm", "w");

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    Image gradient = make_image(image_width, image_height);

    RandomState rng = make_random();
    render_image(&real, &gradient, &rng);
    image_write_bpm(&real, freal);
    image_write_bpm(&gradient, fgradient);
    free_image(&real);
    free_image(&gradient);
}
