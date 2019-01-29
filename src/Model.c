#include <stdlib.h>
#include <math.h>

#include "Model.h"

int getPieceIndex(double v, TypePiecewiseModelParam *param) {
	int a = 0, b = param->size, c;
	if(v<param->startTime[0]) {
		fprintf(stderr, "Error in 'getPieceIndex': value too small\n");
		return 0;
	}
	if(v>param->startTime[param->size]) {
		fprintf(stderr, "Error in 'getPieceIndex': value too high\n");
		return param->size-1;
	}
	while(b > a+1) {
		c = (a+b)/2;
		if(param->startTime[c] < v)
			a = c;
		else
			b = c;
	}
	return a;
}
			

void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param) {
	int i;
	fprintf(f, "size %d\n", param->size);
	for(i=0; i<param->size; i++)
		fprintf(f, "%lf\n%lf %lf %lf\n", param->startTime[i], param->param[i].birth, param->param[i].death, param->param[i].sampl);
	fprintf(f, "%lf\n", param->startTime[param->size]);
}
