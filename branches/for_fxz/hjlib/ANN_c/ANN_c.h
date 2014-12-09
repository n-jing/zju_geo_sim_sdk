#ifndef HJ_ANN_C_INTERFACE_H_
#define HJ_ANN_C_INTERFACE_H_

#ifndef HJ_ANN_C_API
#define HJ_ANN_C_API
#endif

void *HJ_ANN_C_API ANNkd_tree_new(double **pts, int num, int dim);

void HJ_ANN_C_API ANNkd_tree_search(void *ANNkd_tree_handle, double *pt, int num,
									int *idx, double *dist2);

int HJ_ANN_C_API ANNkd_tree_FR_search(void *ANNkd_tree_handle, double *pt, double sq_r, int num,
									  int *idx, double *dist2);

void HJ_ANN_C_API ANNkd_tree_delete(void *ANNkd_tree_handle);

#endif

