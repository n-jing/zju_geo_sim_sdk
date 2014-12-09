#include <ANN/ANN.h>

#ifdef WIN32
#define HJ_ANN_C_API __declspec(dllexport)
#endif

extern "C" {
#include "ANN_c.h"
}

void *ANNkd_tree_new(double **pts, int num, int dim)
{
	return new ANNkd_tree(pts, num, dim);
}

void ANNkd_tree_search(void *ANNkd_tree_handle, double *pt, int num,
					   int *idx, double *dist2)
{
	reinterpret_cast<ANNkd_tree *>(ANNkd_tree_handle)->annkSearch(
		pt, num, idx, dist2);
}

int ANNkd_tree_FR_search(void *ANNkd_tree_handle, double *pt, double sq_r, int num,
					   int *idx, double *dist2)
{
	return reinterpret_cast<ANNkd_tree *>(ANNkd_tree_handle)->annkFRSearch(
		pt, sq_r, num, idx, dist2);
}

void ANNkd_tree_delete(void *ANNkd_tree_handle)
{
	delete reinterpret_cast<ANNkd_tree *>(ANNkd_tree_handle);
}
