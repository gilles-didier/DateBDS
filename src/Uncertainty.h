#ifndef UncertaintyF
#define UncertaintyF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Utils.h"
#include "Tree.h"
#include "Model.h"
#include "Distribution.h"


typedef struct UNCERTAINTY_LR_DATA {
	int nLeaves;
	double logRank;
} TypeUncertaintyLogRankData;


typedef struct UNCERTAINTY_LIKE_PARAMETERS {
    double birth, death, sampl;
    double bmd, bs, bsd;
} TypeUncertaintyLikeParameters;


#ifdef __cplusplus
extern "C" {
#endif
void fillUncertaintyLogRankData(int node, TypeTree *tree, TypeUncertaintyLogRankData *lr);
TypeUncertaintyLikeParameters getUncertaintyLikeParameters(TypeModelParam *param);
TypeDistribution getDistributionUn(int n, double val, TypeTree *tree, TypeModelParam *param);
TypeDistribution getDistribution(int n, int def, TypeTree *tree, TypeModelParam *param);
void fillDistributionAll(int def, TypeTree *tree, TypeModelParam *param, int *todo, TypeDistribution *d);
double getDistributionRoot(double val, TypeTree *tree, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeParameters *param);

#ifdef __cplusplus
}
#endif


#endif
