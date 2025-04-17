#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

#include "hello.h"

float gradient_scaling_factor(long image_width, long image_height) {
    long size = image_width * image_height * 3;
    return 1.f / (float)size;
}

float mse_loss(Image *real, Image *groundtruth) {
    float loss = 0.f;
    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            vec3 radiance_groundtruth;
            vec3 radiance_real;
            image_get(radiance_groundtruth, groundtruth, ir, ic);
            image_get(radiance_real, real, ir, ic);
            vec3 error;
            vec3_sub(error, radiance_real, radiance_groundtruth);
            vec3 squared_error = {powf(error[0], 2), powf(error[1], 2), powf(error[2], 2)};
            float squared_error_avg = squared_error[0] + squared_error[1] + squared_error[2];
            loss += squared_error_avg;
        }
    }
    return loss * gradient_scaling_factor(real->image_width, real->image_height);
}

void mse_loss_deriv(Image *real, Image *groundtruth, GradientImage *gradient, SceneParams *loss_deriv) {
    scene_params_fill(loss_deriv, 0.0);

    SceneParamsPerChannel ppc;
    for (int c = 0; c < 3; c++) {
        ppc.rgb[c] = make_scene_params();
    }

    float factor = gradient_scaling_factor(real->image_width, real->image_height);

    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            vec3 pixel_groundtruth;
            vec3 pixel_real;

            image_get(pixel_groundtruth, groundtruth, ir, ic);
            image_get(pixel_real, real, ir, ic);

            gradient_image_get(&ppc, gradient, ir, ic);

            vec3 error;
            vec3_sub(error, pixel_real, pixel_groundtruth);
            vec3_scale(error, error, 2.f);
            vec3_scale(error, error, factor);

            for (int ch = 0; ch < 3; ch++) {
                scene_params_scale(ppc.rgb[ch], ppc.rgb[ch], error[ch]);
                scene_params_elementwise_add(loss_deriv, loss_deriv, ppc.rgb[ch]);
            }
        }
    }

    for (int ch = 0; ch < 3; ch++) {
        free_scene_params(ppc.rgb[ch]);
    }
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

int main(void) {
    mkdir("temp", 0777); // do nothing if temp already exists

    // FILE *freal = fopen("real.bpm", "w");
    // FILE *fgradient = fopen("gradient.bpm", "w");
    FILE *fgroundtruth = fopen("groundtruth.ppm", "r");

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);
    Image groundtruth = make_image(image_width, image_height);
    image_read_bpm(&groundtruth, fgroundtruth);

    SceneParams *params = make_scene_params();
    params->offset = 0.1f;
    SceneParams *loss_deriv = make_scene_params();

    const int num_epochs = 20;
    for (int epoch = 0; epoch < num_epochs; epoch++) {
        render_image(&real, &gradient, params); // calculate radiance and gradients

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
        image_write_bpm(&real, freal);
    }

    free_scene_params(loss_deriv);
    free_scene_params(params);
    free_image(&real);
    free_gradient_image(&gradient);
    free_image(&groundtruth);
}
