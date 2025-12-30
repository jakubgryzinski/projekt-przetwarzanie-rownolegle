#include <cuda_runtime.h>
#include <stdio.h>

void cpu_sequential(float* input, float* output, int N, int R);
__global__ void gpu_efficient_global(float* input, float* output, int N, int R, int k);
__global__ void gpu_inefficient_global(float* input, float* output, int N, int R, int k);
__global__ void gpu_efficient_shared(float* input, float* output, int N, int R, int k);
__global__ void gpu_inefficient_shared(float* input, float* output, int N, int R, int k);

void run_experiment(int N, int R, int k, int BS);

int main() {
    run_experiment(1024, 2, 1, 16);
    return 0;
}

void run_experiment(int N, int R, int k, int BS) {
    printf("=== Eksperyment: N=%d, R=%d, k=%d, BS=%d ===\n", N, R, k, BS);

    dim3 block(BS, BS);
    dim3 grid(1, 1);
    
    cpu_sequential(NULL, NULL, N, R);
    
    gpu_efficient_global<<<grid, block>>>(NULL, NULL, N, R, k);
    gpu_inefficient_global<<<grid, block>>>(NULL, NULL, N, R, k);
    gpu_efficient_shared<<<grid, block>>>(NULL, NULL, N, R, k);
    gpu_inefficient_shared<<<grid, block>>>(NULL, NULL, N, R, k);
    
    cudaDeviceSynchronize();
    printf("\n");
}

void cpu_sequential(float* input, float* output, int N, int R) {
    printf("CPU sekwencyjny: N=%d, R=%d\n", N, R);
}

__global__ void gpu_efficient_global(float* input, float* output, int N, int R, int k) {
    printf("GPU efektywny global: N=%d, R=%d, k=%d\n", N, R, k);
}

__global__ void gpu_inefficient_global(float* input, float* output, int N, int R, int k) {
    printf("GPU nieefektywny global: N=%d, R=%d, k=%d\n", N, R, k);
}

__global__ void gpu_efficient_shared(float* input, float* output, int N, int R, int k) {
    printf("GPU efektywny shared: N=%d, R=%d, k=%d\n", N, R, k);
}

__global__ void gpu_inefficient_shared(float* input, float* output, int N, int R, int k) {
    printf("GPU nieefektywny shared: N=%d, R=%d, k=%d\n", N, R, k);
}