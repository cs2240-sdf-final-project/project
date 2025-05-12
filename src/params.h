#ifndef PARAMS_H
#define PARAMS_H

#include <cassert>
#include "linmath.h"

extern int number_of_scene_params;

typedef struct SceneParams SceneParams;
typedef struct SceneParamsPerChannel SceneParamsPerChannel;

extern SceneParams *uninit_scene_params(void);
extern void free_scene_params(SceneParams *params);

struct SceneParamsPerChannel {
    SceneParams *rgb[3];
};

static inline const float *float_pointer_from_params(const SceneParams *out) {
    return (float *)out;
}

static inline float *float_pointer_from_params(SceneParams *out) {
    return (float *)out;
}

static inline float scene_parameter_get(const SceneParams *params, long p) {
    const float *raw_params = float_pointer_from_params(params);
    return raw_params[p];
}

static inline void scene_params_set(SceneParams *params, long p, float value) {
    float *raw_params = float_pointer_from_params(params);
    raw_params[p] = value;
}

static inline void params_nan_to_num(SceneParamsPerChannel *ppc, float num) {
    for (int ch = 0; ch < 3; ch++) {
        for (long p = 0; p < number_of_scene_params; p++) {
            float param = scene_parameter_get(ppc->rgb[ch], p);
            if (isnan(param)) {
                scene_params_set(ppc->rgb[ch], p, num);
            }
        }
    }
}

static inline SceneParams *index_channel(SceneParamsPerChannel *ppc, long channel) {
    assert(channel < 3);
    return ppc->rgb[channel];
}

static inline void params_per_channel_add_assign(SceneParamsPerChannel *ppc, const SceneParamsPerChannel *to_add) {
    for (int c = 0; c < 3; c++) {
        float *apc = float_pointer_from_params(ppc->rgb[c]);
        const float *bpc = float_pointer_from_params(to_add->rgb[c]);
        for (int i = 0; i < number_of_scene_params; i++) {
            apc[i] += bpc[i];
        }
    }
}

static inline float *ppc_get_channel(SceneParamsPerChannel *ppc, long ch) {
    assert(ch < 3);
    return float_pointer_from_params(ppc->rgb[ch]);
}

static inline const float *ppc_get_channel(const SceneParamsPerChannel *ppc, long ch) {
    assert(ch < 3);
    return float_pointer_from_params(ppc->rgb[ch]);
}

static inline void scene_params_copy(SceneParams *out, const SceneParams *params) {
    float *raw_out = float_pointer_from_params(out);
    const float *raw_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out[i] = raw_params[i];
    }
}

static inline void scene_params_elementwise_add(SceneParams *out_params, const SceneParams *a, const SceneParams *b) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    const float *raw_b = float_pointer_from_params(b);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] + raw_b[i];
    }
}

static inline void scene_params_elementwise_mul(SceneParams *out_params, const SceneParams *a, const SceneParams *b) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    const float *raw_b = float_pointer_from_params(b);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] * raw_b[i];
    }
}

static inline void outer_product_add_assign(SceneParamsPerChannel *ppc, const SceneParams *params, const vec3 rgb) {
    const float *raw_params = float_pointer_from_params(params);
    for (int c = 0; c < 3; c++) {
        float *pc = float_pointer_from_params(ppc->rgb[c]);
        for (int i = 0; i < number_of_scene_params; i++) {
            pc[i] += raw_params[i] * rgb[c];
        }
    }
}

static inline void scene_params_scale(SceneParams *out_params, const SceneParams *a, float scale_by) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] * scale_by;
    }
}

static inline void scene_params_fill(SceneParams *params, float fill_with) {
    float *raw_out_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = fill_with;
    }
}

static inline SceneParamsPerChannel *make_ppc(void) {
    SceneParamsPerChannel *ppc = new SceneParamsPerChannel;
    for (int c = 0; c < 3; c++) {
        ppc->rgb[c] = uninit_scene_params();
    }
    return ppc;
}

static inline void set_ppc(SceneParamsPerChannel *ppc, long param, const vec3 value) {
    for (int ch = 0; ch < 3; ch++) {
        scene_params_set(ppc->rgb[ch], param, value[ch]);
    }
}

static inline void ppc_scale(SceneParamsPerChannel *ppc, vec3 scale_by) {
    for (int ch = 0; ch < 3; ch++) {
        float *floats = ppc_get_channel(ppc, ch);
        for (int p = 0; p < number_of_scene_params; p++) {
            floats[p] *= scale_by[ch];
        }
    }
}

static inline void fill_ppc(SceneParamsPerChannel *ppc, float fill_with) {
    for (int ch = 0; ch < 3; ch++) {
        scene_params_fill(ppc->rgb[ch], fill_with);
    }
}

static inline void free_ppc(SceneParamsPerChannel *ppc) {
    for (int ch = 0; ch < 3; ch++) {
        free_scene_params(ppc->rgb[ch]);
    }
    delete ppc;
}

static inline SceneParams *params_from_float_pointer(const float *params) {
    return (SceneParams *)params;
}

#endif // PARAMS_H
