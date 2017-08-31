/* Computes a step in the integration changing both p and r */

__global__ void step_type1(long n, double a1, double a2, double *r_gpu, double *p_gpu, double *f_gpu)
{

    long tid;

    tid=threadIdx.x+blockIdx.x*blockDim.x;

    while (tid<n)
    {
        p_gpu[tid]=p_gpu[tid]+a1*f_gpu[tid];
        r_gpu[tid]=r_gpu[tid]+a2*p_gpu[tid];

        tid+=blockDim.x*gridDim.x;
    }

    return;

}

