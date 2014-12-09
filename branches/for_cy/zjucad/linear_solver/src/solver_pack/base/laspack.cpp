#include <string.h>

#include <iostream>

#include "laspack.h"

using namespace std;

struct laspack_iter_solver
{
	const char *name_;
	const IterProcType proc_;
};

struct laspack_precond
{
	const char *name_;
	const PrecondProcType proc_;
};

const laspack_iter_solver laspack_iter_solvers[] = {
	{"Jacobi", JacobiIter},
	{"SORForw", SORForwIter},
	{"SORBackw", SORBackwIter},
	{"SSOR", SSORIter},
	{"Chebyshev", ChebyshevIter},
	{"CG", CGIter},
	{"CGN", CGNIter},
	{"GMRES", GMRESIter},
	{"BiCG", BiCGIter},
	{"QMR", QMRIter},
	{"CGS", CGSIter},
	{"BiCGSTAB", BiCGSTABIter}
};

const laspack_precond laspack_preconds[] = {
	{"Jacobi", JacobiPrecond},
	{"SSOR", SSORPrecond},
	{"ILU", ILUPrecond}
};

size_t get_iter_type(const char *name[])
{
	if(!name)
		return sizeof(laspack_iter_solvers)/sizeof(laspack_iter_solver);
	for(size_t i = 0; i < get_iter_type(0); ++i)
		name[i] = laspack_iter_solvers[i].name_;
	return 0;
}

IterProcType get_iter_proc(const char *name)
{
	for(size_t i = 0; i < sizeof(laspack_iter_solvers)/sizeof(laspack_iter_solver); ++i) {
		if(!strcmp(laspack_iter_solvers[i].name_, name))
			return laspack_iter_solvers[i].proc_;
	}
	return 0;
}

size_t get_precond_type(const char *name[])
{
	if(!name)
		return sizeof(laspack_preconds)/sizeof(laspack_precond);
	for(size_t i = 0; i < get_precond_type(0); ++i)
		name[i] = laspack_preconds[i].name_;
	return 0;
}

PrecondProcType get_precond_proc(const char *name)
{
	for(size_t i = 0; i < sizeof(laspack_preconds)/sizeof(laspack_precond); ++i) {
		if(!strcmp(laspack_preconds[i].name_, name))
			return laspack_preconds[i].proc_;
	}
	return 0;
}
