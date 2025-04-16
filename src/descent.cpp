#include <stdio.h>
#include <stdlib.h>
#include "hello.h"
#include <sys/stat.h>

float mse_loss(Image *real, Image *groundtruth) {
    float loss = 0.f;
    float count = 0.f;
    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            vec3 radiance_groundtruth;
            vec3 radiance_real;
            image_get(radiance_groundtruth, groundtruth, ir, ic);
            image_get(radiance_real, real, ir, ic);
            vec3 error;
            vec3_sub(error, radiance_real, radiance_groundtruth);
            vec3 squared_error = {powf(error[0], 2), powf(error[1], 2), powf(error[2], 2)};
            float squared_error_avg = (1.f / 3.f) * (squared_error[0] + squared_error[1] + squared_error[2]);
            loss += squared_error_avg;
            count += 1.f;
        }
    }
    loss /= count;
    return loss;
}

void mse_loss_deriv(Image *real, Image *groundtruth, Image *gradient, SceneParams *loss_deriv) {
    float deriv = 0.f;
    float count = 0.f;
    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            vec3 pixel_groundtruth;
            vec3 pixel_real;
            vec3 pixel_gradient;

            image_get(pixel_groundtruth, groundtruth, ir, ic);
            image_get(pixel_real, real, ir, ic);
            image_get(pixel_gradient, gradient, ir, ic);

            vec3 error;
            vec3_sub(error, pixel_real, pixel_groundtruth);
            vec3 twice_error;
            vec3_scale(twice_error, error, 2.f);
            vec3 pixel_deriv = {twice_error[0] * pixel_gradient[0],
                twice_error[1] * pixel_gradient[1],
                twice_error[2] * pixel_gradient[2]};

            float pixel_deriv_avg = (1.f / 3.f) * (pixel_deriv[0] + pixel_deriv[1] + pixel_deriv[2]);
            deriv += pixel_deriv_avg;
            count += 1.f;
        }
    }
    deriv /= count;
    return deriv;
}

void gradient_step(SceneParams *params, const SceneParams *deriv, float learning_rate = 500.f) {
    // TODO: implement a better optimizer like Adam or something
    // *param = *param - learning_rate * deriv;
    SceneParams *scratch = make_scene_params();
    scene_params_fill(scratch, -learning_rate);
    scene_params_elementwise_mul(scratch, scratch, deriv);
    scene_params_elementwise_add(params, params, scratch);
    free_scene_params(scratch);
}

int main(int argc, char *argv[]) {
    mkdir("temp", 0777); // do nothing if temp already exists

    // FILE *freal = fopen("real.bpm", "w");
    // FILE *fgradient = fopen("gradient.bpm", "w");
    FILE *fgroundtruth = fopen("groundtruth.ppm", "r");

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    Image gradient = make_image(image_width, image_height);
    Image groundtruth = make_image(image_width, image_height);
    image_read_bpm(&groundtruth, fgroundtruth);

    SceneParams *params = make_scene_params();
    params->offset = 0.1f;
    SceneParams *loss_deriv = make_scene_params();

    const int num_epochs = 20;
    for (int epoch = 0; epoch < num_epochs; epoch++) {
        RandomState rng = make_random();
        render_image(&real, &gradient, &rng, &params); // calculate radiance and gradients

        // Compute loss and derivative of loss
        float loss = mse_loss(&real, &groundtruth);
        mse_loss_deriv(&real, &groundtruth, &gradient, loss_deriv);

        printf("loss: %f\n", loss);
        // printf("deriv: %f\n", loss_deriv); // TODO: print this struct

        // Gradient step
        gradient_step(params, loss_deriv, 500.0f);

        // printf("param: %f\n", param);

        // Write each frame to a file
        char fn_real[256];
        snprintf(fn_real, sizeof(fn_real), "temp/real_%04d.ppm", epoch);
        FILE *freal = fopen(fn_real, "w");
        char fn_gradient[256];
        snprintf(fn_gradient, sizeof(fn_gradient), "temp/gradient_%04d.ppm", epoch);
        FILE *fgradient = fopen(fn_gradient, "w");
        image_write_bpm(&real, freal);
        image_write_bpm(&gradient, fgradient);
    }

    free_scene_params(loss_deriv);
    free_scene_params(params);
    free_image(&real);
    free_image(&gradient);
    free_image(&groundtruth);
}
