#include "common.h"

#include <fstream>

using namespace std;

void dump(size_t iter_count, size_t dump_step, const char *dump_pref,
		  const double *x, int32_t nrow, int32_t ncol)
{
	if(dump_step > 0 && (iter_count%dump_step == 0)) {
		char name[256];
		sprintf(name, "%s-%05ld.x.mat", dump_pref, iter_count);
		ofstream ofs(name, ofstream::binary);
		ofs.write((const char *)&nrow, sizeof(int));
		ofs.write((const char *)&ncol, sizeof(int));
		ofs.write((const char *)x, sizeof(double)*nrow*ncol);
	}
}
