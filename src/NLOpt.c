#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <nlopt.h>

#include "Utils.h"

#include "NLOpt.h"


//#define NLOPT_ALGO NLOPT_GN_ISRES
//#define NLOPT_ALGO NLOPT_GN_ESCH
#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_AUGLAG


#define MINVAL 0.01
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define TOLERANCE_CONSTRAINT 0.000000001
#define TOLERANCE_OPTIM 0.001





#define TAG_TRI "TRI"
#define TAG_TOL "TOL"
#define TAG_ITE "ITE"
#define SIZE_TAG 20
#define SIZE_VAL 100

void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    fprintf(f, ":%s %d\n", TAG_TRI, option->trials);
    fprintf(f, ":%s %lE\n", TAG_TOL, option->tolOptim);
    fprintf(f, ":%s %d\n", TAG_ITE, option->maxIter);
}

void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    char c, tag[SIZE_TAG+1], val[SIZE_VAL+1];
    for(c=fgetc(f); c!=EOF && isspace(c); c=fgetc(f));
    while(c == ':') {
        int i;
        c=fgetc(f);
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_TAG; c=fgetc(f))
            tag[i++] = c;
        tag[i] = '\0';
        if(i>=SIZE_TAG) {
            fprintf(stderr, "Error when reading an optimizer options file - Tag too long:\n%s...\n", tag);
            exit(1);
        }
        for(; c!=EOF && isspace(c); c=fgetc(f));
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_VAL; c=fgetc(f))
            val[i++] = c;
        val[i] = '\0';
        if(i>=SIZE_VAL) {
            fprintf(stderr, "Error when reading an optimizer options file - value too long:\n%s...\n", val);
            exit(1);
        }
        if(strcmp(tag, TAG_TRI) == 0)
            option->trials = atoi(val);
        if(strcmp(tag, TAG_TOL) == 0)
            option->tolOptim = atof(val);
        if(strcmp(tag, TAG_ITE) == 0)
            option->maxIter = atoi(val);
        for(; c!=EOF && isspace(c); c=fgetc(f));
    }
}

void fprintNLoptOption(FILE *f, TypeNLOptOption *option) {
    fprintf(f, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void sprintNLoptOption(char *buffer, TypeNLOptOption *option) {
    buffer += sprintf(buffer, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}


int maximizeNLOPT(unsigned size, double *estim, double *lowerBounds, double *upperBounds, double *max, TypeFuncNLOPT f, void *data, TypeNLOptOption *option) {
    double *x, maxTmp;
    nlopt_opt opt;
    int result, t;
    x = (double*) malloc(size*sizeof(double));
    opt = nlopt_create(NLOPT_ALGO, size); /* algorithm and dimensionality */
//	nlopt_add_inequality_constraint(opt, ratesConstraint, NULL, 1e-8);
    nlopt_set_lower_bounds(opt, lowerBounds);
    nlopt_set_upper_bounds(opt, upperBounds);
    nlopt_set_max_objective(opt, f, data);
    nlopt_set_xtol_abs1(opt, option->tolOptim);
    nlopt_set_maxeval(opt, option->maxIter);
    *max = -INFTY;
    for(t=0; t<option->trials; t++) {
		unsigned i;
		for(i=0; i<size; i++)
			x[i] = lowerBounds[i]+UNIF_RAND*(upperBounds[i]-lowerBounds[i]);
        if(((result = nlopt_optimize(opt, x, &maxTmp)) >= 0) && maxTmp > *max) {
			for(i=0; i<size; i++)
				estim[i] = x[i];
			*max = maxTmp;
        }
    }
    free((void*)x);
    nlopt_destroy(opt);
    return result;
}
