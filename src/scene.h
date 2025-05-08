#ifndef SCENE_H
#define SCENE_H

typedef struct {
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float emissive[3];
    bool isReflected;
    float shininess;
} SdfResult;

typedef struct SceneContext SceneContext;

SceneContext *make_scene_context();
void free_scene_context(SceneContext *ctx);

typedef struct SceneParams SceneParams;

void scene_params_init(SceneParams *params, const SceneContext *ctx);
SceneParams *uninit_scene_params();
void free_scene_params(SceneParams *params);

float scene_consistency_loss(const SceneParams *params);
void scene_consistency_gradient(const SceneParams *params, SceneParams *gradient_out);
float directional_derivative(const vec3 origin, const vec3 direction, float t, const SceneParams *params, const SceneContext *ctx);
void scene_sample(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *res);
float get_normal_from(vec3 out_normal, const vec3 pos, const SceneParams *params, const SceneContext *ctx);
float scene_sample_sdf(const vec3 pos, const SceneParams *params, const SceneContext *ctx);

#endif // SCENE_H
