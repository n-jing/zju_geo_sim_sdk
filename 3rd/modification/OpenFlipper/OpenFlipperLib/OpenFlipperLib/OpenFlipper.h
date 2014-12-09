#ifndef HJ_OPENFLIPPER_LIB_H_
#define HJ_OPENFLIPPER_LIB_H_

#ifdef OpenFlipperLib_EXPORTS
  #ifdef _MSC_VER
    #define OpenFlipperLibAPI __declspec(dllexport)
  #else
    #define OpenFlipperLibAPI __attribute__((dllexport))
  #endif
#else
  #define OpenFlipperLibAPI
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct tri_mesh
{
	double *vertex;
	int *face;
	int vert_num, face_num;
};

OpenFlipperLibAPI
int IsotropicRemesher0(
	const tri_mesh *input,
	double edge_length,
	int *vert_num, int *face_num,
	void **ctx);

OpenFlipperLibAPI
int IsotropicRemesher1(
	void *ctx,
	tri_mesh *output);

#ifdef __cplusplus
}
#endif

#endif
