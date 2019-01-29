#ifndef NLOptF
#define NLOptF

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*TypeFuncNLOPT)(unsigned, const double *, double *, void *);

typedef struct NLOPT_OPTION {
    double tolOptim;
    int trials, maxIter;
} TypeNLOptOption;


int maximizeNLOPT(unsigned size, double *estim, double *lowerBounds, double *upperBounds, double *max, TypeFuncNLOPT f, void *data, TypeNLOptOption *option);

void fprintNLoptOption(FILE *f, TypeNLOptOption *option);
void sprintNLoptOption(char *buffer, TypeNLOptOption *option);
void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option);
void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option);
#ifdef __cplusplus
}
#endif

#endif
