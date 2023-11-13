#include <cuda_runtime.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <thread>

__global__ void add(float* x, float* y, float* z, int n)
{
    // 获取全局索引
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // 步长
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
    {
        z[i] = x[i] * y[i];
    }
}

/*
int main(void)
{
    int N = 1e8;
    int nBytes = N * sizeof(float);

    // 申请托管内存
    float* x, * y, * z1, *z2;
    cudaMallocManaged((void**)&x, nBytes);
    cudaMallocManaged((void**)&y, nBytes);
    cudaMallocManaged((void**)&z1, nBytes);

    // 初始化数据
    for (int i = 0; i < N; ++i)
    {
        x[i] = 11110.0;
        y[i] = 22220.0;
    }

    clock_t start1 = clock();
    for (int i = 0; i < N; i++) {
        z2[i] = x[i] * y[i];
    }
    clock_t end1 = clock();
    clock_t start2 = clock();
    // 定义kernel的执行配置
    dim3 blockSize(256);
    dim3 gridSize((N + blockSize.x - 1) / blockSize.x);
    // 执行kernel
    add << < gridSize, blockSize >> > (x, y, z1, N);
    // 同步device 保证结果能正确访问
    cudaDeviceSynchronize();
    // 检查执行结果
    clock_t end2 = clock();
    
    // output
    double t1 = ((double)(end1 - start1)) / CLOCKS_PER_SEC;
    double t2 = ((double)(end2 - start2)) / CLOCKS_PER_SEC;
    std::cout << "time (CPU): " << t1 * 1000 << " ms" << std::endl;
    std::cout << "time (GPU): " << t2 * 1000 << " ms" << std::endl;

    // 释放内存
    cudaFree(x);
    cudaFree(y);
    cudaFree(z1);
    cudaFree(z2);

    return 0;
}
*/