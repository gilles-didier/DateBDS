#ifndef ModelF
#define ModelF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Tree.h"

#define RINFTY 1E99

typedef struct MODEL_PARAM {
	double birth, death, sampl;
} TypeModelParam;

typedef struct PIECEWISE_MODEL_PARAM {
	int size;
	double *startTime;
	TypeModelParam *param;
} TypePiecewiseModelParam;

#ifdef __cplusplus
extern "C" {
#endif

int getPieceIndex(double v, TypePiecewiseModelParam *param);
void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param);

#ifdef __cplusplus
}
#endif

#endif
