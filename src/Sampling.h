#ifndef SamplingF
#define SamplingF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "Tree.h"
#include "Model.h"

#ifdef __cplusplus
extern "C" {
#endif

void sampleTimes(double eps, TypeTree *tree, TypeModelParam *param, gsl_rng *rg);

#ifdef __cplusplus
}
#endif



#endif
