#ifdef __cplusplus
extern "C" {
#endif

void BetaPDF(int nx, double *x, double mean, double var, double *pdf,
                     int *pdfBound);
double scaledBetaDistribution(double x, double dx, double alpha, double beta);
double lnGamma(double xx);

#ifdef __cplusplus
}
#endif
