#ifndef DistributionF
#define DistributionF

#include <stdlib.h>
#include <stdio.h>

typedef struct DENSITY_ITEM {
    double val, dens;
} TypeDistributionItem;

typedef struct DENSITY {
    int size;
    TypeDistributionItem *item;
} TypeDistribution;

#ifdef __cplusplus
extern "C" {
#endif


int compareDistributionItem(const void* a, const void* b);
TypeDistribution meanDistribution(TypeDistribution d1, double w1, TypeDistribution d2, double w2);
void fprintDistribution(FILE *f, TypeDistribution d);
double sumDistribution(TypeDistribution d);
void fprintDerive(FILE *f, TypeDistribution d);
double getMean(TypeDistribution d);
double getMedian(TypeDistribution d);
double getMeanDens(TypeDistribution d);
double getMedianDens(TypeDistribution d);
double getQuantileInf(TypeDistribution d, double q, double x);
double getQuantileSup(TypeDistribution d, double q, double x);
void deriveDistribution(TypeDistribution *d);

#ifdef __cplusplus
}
#endif



#endif
