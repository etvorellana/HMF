#ifndef openmpLibs_h
#define openmpLibs_h
void WaterBag(long, long *, double, double, double *, double *);
void KineticEnergy(long, double *, double *);
void PotentialEnergy(long, double *, double *, double, double);
void Force(long, double *, double *, double *, double *);
void Integration(long, double, double *, double *, double *, double *, double *);
float ran2(long *);
#endif
