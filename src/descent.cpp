#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

#include "hello.h"
#include "sim_random.h"

float gradient_scaling_factor(long image_width, long image_height) {
    long size = image_width * image_height * 3;
    return 1.f / (float)size;
}

float mse_loss(const Image *real, const Image *groundtruth) {
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

void mse_loss_deriv(const Image *real, const Image *groundtruth, GradientImage *gradient, SceneParams *loss_deriv) {
    scene_params_fill(loss_deriv, 0.0);

    SceneParamsPerChannel ppc;
    for (int c = 0; c < 3; c++) {
        ppc.rgb[c] = uninit_scene_params();
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
            vec3_scale(error, error, factor * 2.0f);

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

float total_loss(const Image *real, const Image *groundtruth, const SceneParams *params) {
    float loss_image = mse_loss(real, groundtruth);
    float loss_scene_consistency = scene_consistency_loss(params);
    return loss_image + loss_scene_consistency;
}

void total_loss_deriv(const Image *real, Image *groundtruth, GradientImage *gradient, SceneParams *loss_deriv_out, SceneParams *scratch, const SceneParams *params) {
    scene_consistency_gradient(params, scratch);
    mse_loss_deriv(real, groundtruth, gradient, loss_deriv_out);
    scene_params_elementwise_add(loss_deriv_out, loss_deriv_out, scratch);
}

void gradient_step(SceneParams *params, const SceneParams *deriv, float learning_rate, SceneParams *scratch) {
    // TODO: implement a better optimizer like Adam or something
    // *param = *param - learning_rate * deriv;
    scene_params_fill(scratch, -learning_rate);
    scene_params_elementwise_mul(scratch, scratch, deriv);
    scene_params_elementwise_add(params, params, scratch);
}

int main(void) {
    mkdir("descent-sequence", 0777); // do nothing if temp already exists

    FILE *fgroundtruth = fopen("groundtruth.ppm", "r");
    assert(fgroundtruth);

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    GradientImage gradient = make_gradient_image(image_width, image_height);
    Image groundtruth = make_image(image_width, image_height);
    image_read_bpm(&groundtruth, fgroundtruth);

    SceneContext *ctx = make_scene_context();

    RandomState *rng = make_random();
    SceneParams *params = uninit_scene_params();
    scene_params_init(params, ctx);
    SceneParams *loss_deriv = uninit_scene_params();
    SceneParams *scratch = uninit_scene_params();

    const float learning_rate = 1e-1f;

    const int num_epochs = 1000;
    for (int epoch = 0; epoch < num_epochs; epoch++) {
        render_image_phong(&real, &gradient, params,ctx, rng); // calculate radiance and gradients

        // Compute loss and derivative of loss
        float loss = total_loss(&real, &groundtruth, params);
        total_loss_deriv(&real, &groundtruth, &gradient, loss_deriv, scratch, params);

        printf("loss: %f\n", loss);
        fflush(stdout);
        printf("deriv"); // TODO: print this struct
        for (int p = 0; p < number_of_scene_params; p++) {
            printf("%f ", scene_parameter_get(loss_deriv, p)); // TODO: print this struct
        }
        printf("\n"); // TODO: print this struct
        fflush(stdout);

        // Gradient step
        gradient_step(params, loss_deriv, learning_rate, scratch);

        // Write each frame to a file
        char fn_real[256];
        snprintf(fn_real, sizeof(fn_real), "descent-sequence/real_%04d.ppm", epoch);
        FILE *freal = fopen(fn_real, "w");
        image_write_ppm(&real, freal);
        fclose(freal);
    }

    free_scene_context(ctx);
    free_random(rng);
    free_scene_params(scratch);
    free_scene_params(loss_deriv);
    free_scene_params(params);
    free_image(&real);
    free_gradient_image(&gradient);
    free_image(&groundtruth);
    fclose(fgroundtruth);
}
