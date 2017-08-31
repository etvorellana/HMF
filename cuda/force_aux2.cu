
__global__ void force_aux2(long n, double mx, double my, double *r_gpu, double *f_gpu, float *sinr_gpu, float *cosr_gpu)
{
    long tid;

    tid=threadIdx.x+blockIdx.x*blockDim.x;

    while (tid<n)
    {
        f_gpu[tid]=cosr_gpu[tid]*my-sinr_gpu[tid]*mx;
        tid+=blockDim.x*gridDim.x;
    };

    return;

}
