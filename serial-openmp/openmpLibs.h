
#define	pi	3.14159265359
#define dpi	6.28318530718
#define B0   0.675603595979828813
#define B1  -0.175603595979828813
#define D0   1.35120719195965763
#define D1  -1.70241438391931525

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define RMAX 1.0E-5

void WaterBag(long, long *, double, double, double *, double *);
void KineticEnergy(long, double *, double *);
void PotentialEnergy(long, double *, double *, double, double);
void Force(long, double *, double *, double *, double *);
void Integration(long, double, double *, double *, double *, double *, double *);
float ran2(long *idum);
