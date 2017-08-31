/* Computes a step in the integration changing only p */

__global__ void step_type2(long n, double a, double *r_gpu, double *p_gpu, double *f_gpu)
{

    long tid;

    tid=threadIdx.x+blockIdx.x*blockDim.x;

    while (tid<n)
    {
        p_gpu[tid]=p_gpu[tid]+a*f_gpu[tid];

        tid+=blockDim.x*gridDim.x;
    }

    return;

}

