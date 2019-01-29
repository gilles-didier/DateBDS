#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "Distribution.h"



int compareDistributionItem(const void* a, const void* b) {
    if(((TypeDistributionItem*)a)->val>((TypeDistributionItem*)b)->val)
        return 1;
    if(((TypeDistributionItem*)a)->val<((TypeDistributionItem*)b)->val)
        return -1;
    return 0;
}

TypeDistribution meanDistribution(TypeDistribution d1, double w1, TypeDistribution d2, double w2) {
	TypeDistribution d;
	if(d1.size != d2.size) {
		d.size = 0;
		d.item = NULL;
		return d;
	}
	return d;
}

double sumDistribution(TypeDistribution d) {
	int i;
	double sum;
	if(d.size == 0)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	sum = (d.item[1].val-d.item[0].val)*d.item[0].dens;
	for(i=1; i<d.size-1; i++)
		sum += ((d.item[i].val-d.item[i-1].val)/2.+(d.item[i+1].val-d.item[i].val)/2.)*d.item[i].dens;
	sum += (d.item[d.size-1].val-d.item[d.size-2].val)*d.item[d.size-1].dens;
	return sum;
}

void fprintDistribution(FILE *f, TypeDistribution d) {
	int i;
	for(i=0; i<d.size; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, d.item[i].dens);
}

void fprintDerive(FILE *f, TypeDistribution d) {
	int i;
	for(i=0; i<d.size-1 && d.item[i].dens == 0.; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, 0.);
	for(; i<d.size-1 && d.item[i].dens > 0.; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, (d.item[i+1].dens-d.item[i].dens)/(d.item[i+1].val-d.item[i].val));
	for(; i<d.size-1; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, 0.);
	fprintf(f, "%lf\t%le\n", d.item[d.size-1].val, (d.item[d.size-1].dens-d.item[d.size-2].dens)/(d.item[d.size-1].val-d.item[d.size-2].val));
	fprintf(f, "%lf\t%le\n", d.item[d.size-1].val, 0.);
}

//void deriveDistribution(TypeDistribution *d) {
	//if(d->size>0) {
		//int i;
		//for(i=0; i<d->size-1; i++)
			//d->item[i].dens = (d->item[i+1].dens-d->item[i].dens)/(d->item[i+1].val-d->item[i].val);
		//d->item[d->size-1].dens = (d->item[d->size-1].dens-d->item[d->size-2].dens)/(d->item[d->size-1].val-d->item[d->size-2].val);
////		d->item[d->size-1].dens = 0.;
////		d->item[0].dens = 0.;
	//}
//}

void deriveDistribution(TypeDistribution *d) {
//printf("dist %d\n", d->size);
	if(d->size>0) {
		int i;
		double scale = (d->item[d->size-1].val-d->item[0].val)/(d->item[d->size-2].val-d->item[0].val);
//printf("Mx %.2le %.2le\n", d->item[0].val+(d->item[d->size-2].val-d->item[0].val)*scale, d->item[d->size-1].val);
		for(i=0; i<d->size-1; i++) {
			d->item[i].dens = (d->item[i+1].dens-d->item[i].dens)/(d->item[i+1].val-d->item[i].val);
			d->item[i].val = d->item[0].val+(d->item[i].val-d->item[0].val)*scale;
//printf("%d (%.2le, %.2le)\n", i, d->item[i].val, d->item[i].dens);
		}
//		d->item[d->size-1].dens = (d->item[d->size-1].dens-d->item[d->size-2].dens)/(d->item[d->size-1].val-d->item[d->size-2].val);
		d->item[d->size-1].dens = 0.;
//		d->item[0].dens = 0.;
	}
}
/*
void deriveDistribution(TypeDistribution *d) {
	if(d->size>1) {
		int i;
		for(i=1; i<d->size-1; i++)
//			d->item[i].dens = (d->item[i].dens-d->item[i-1].dens)/(d->item[i].val-d->item[i-1].val);
//			d->item[i].dens = (d->item[i].dens-d->item[i-1].dens);
			d->item[i].dens = (d->item[i+1].dens-d->item[i].dens);
		d->item[0].dens = d->item[1].dens;
	}
}
*/

double getMean(TypeDistribution d) {
	int i;
	double sum=0.;
	for(i=0; i<d.size-1; i++)
		sum += d.item[i].val*(d.item[i+1].dens- d.item[i].dens);
	return sum;
//	return sumDistribution(d)+d.item[d.size-1].val+d.item[0].val;
}

double getMedianDens(TypeDistribution d) {
	int i;
	for(i=0; i<d.size && d.item[i].dens <= 0.5; i++)
		;
	return d.item[i-1].val;
}

double getMeanDens(TypeDistribution d) {
	int i;
	double sum=0.;
	for(i=0; i<d.size; i++)
		sum += d.item[i].val*d.item[i].dens;
	return sum;
//	return sumDistribution(d)+d.item[d.size-1].val+d.item[0].val;
}

double getMedian(TypeDistribution d) {
	int i;
	for(i=0; i<d.size && d.item[i].dens <= 0.5; i++)
		;
	return d.item[i-1].val;
}


double getQuantileInf(TypeDistribution d, double q, double x) {
	int i, ix;
	double sum;
	if(d.size == 0)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	for(ix=0; ix<d.size-1 && x<d.item[ix+1].val; ix++)
		;
	for(i=ix; i>0 && sum <= q; i--)
		sum += ((d.item[i].val-d.item[i-1].val)/2.+(d.item[i+1].val-d.item[i].val)/2.)*d.item[i].dens;
	return d.item[i-1].val;
}

double getQuantileSup(TypeDistribution d, double q, double x) {
	int i, ix;
	double sum;
	if(d.size == 0)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	for(ix=1; ix<d.size && x>d.item[ix-1].val; ix++)
		;
	for(i=ix; i<d.size && sum <= q; i++)
		sum += ((d.item[i].val-d.item[i-1].val)/2.+(d.item[i+1].val-d.item[i].val)/2.)*d.item[i].dens;
	return d.item[i-1].val;
}	
