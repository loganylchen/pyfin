#include "dtw_c_api.h"
#include "dtw.hpp"
#include "cuda_utils.hpp"
#include <cuda_runtime.h>
#include <cstdio>

// 复用 OpenDBA 的核函数启动逻辑，仅适配单对序列计算
int opendba_dtw_cuda(
    const float *seq1, size_t len1,
    const float *seq2, size_t len2,
    int use_open_start,
    int use_open_end,
    float *out_distance)
{
    // 1. 输入校验
    if (!seq1 || !seq2 || !out_distance || len1 == 0 || len2 == 0)
    {
        fprintf(stderr, "Invalid input parameters\n");
        return -1;
    }

    // 2. 设备内存分配（复用 OpenDBA 的内存对齐逻辑）
    float *d_seq1, *d_seq2, *d_dtw_cost, *d_new_dtw_cost;
    unsigned char *d_path_matrix;
    float *d_pairwise_dist;
    size_t path_mem_pitch;

    // 分配序列内存
    CUDA_CHECK(cudaMalloc(&d_seq1, len1 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_seq2, len2 * sizeof(float)));
    // 分配 DTW 计算所需临时内存（参考 OpenDBA 原版逻辑）
    CUDA_CHECK(cudaMalloc(&d_dtw_cost, len2 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_new_dtw_cost, len2 * sizeof(float)));
    CUDA_CHECK(cudaMallocPitch(&d_path_matrix, &path_mem_pitch, len2 * sizeof(unsigned char), len1));
    CUDA_CHECK(cudaMalloc(&d_pairwise_dist, sizeof(float)));

    // 3. 主机→设备数据拷贝
    CUDA_CHECK(cudaMemcpy(d_seq1, seq1, len1 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_seq2, seq2, len2 * sizeof(float), cudaMemcpyHostToDevice));
    // 初始化临时内存
    CUDA_CHECK(cudaMemset(d_dtw_cost, 0, len2 * sizeof(float)));
    CUDA_CHECK(cudaMemset(d_new_dtw_cost, 0, len2 * sizeof(float)));
    CUDA_CHECK(cudaMemset(d_path_matrix, 0, path_mem_pitch * len1));
    CUDA_CHECK(cudaMemset(d_pairwise_dist, 0, sizeof(float)));

    // 4. 启动 OpenDBA 原版 DTW 核函数（参数严格对齐）
    int max_threads = getMaxThreadsPerDevice(0); // 复用 OpenDBA 的线程数获取函数
    dim3 thread_block(max_threads, 1, 1);
    size_t shared_mem = thread_block.x * 3 * sizeof(float); // 复用 OpenDBA 的共享内存计算

    DTWDistance<<<1, thread_block, shared_mem>>>(
        d_seq1, len1,
        d_seq2, len2,
        0, 0,          // 单序列对，index/offset 设为 0
        nullptr, 0, 0, // 多序列相关参数置空/0
        nullptr,
        d_dtw_cost,
        d_new_dtw_cost,
        d_path_matrix,
        path_mem_pitch,
        d_pairwise_dist,
        use_open_start,
        use_open_end);
    CUDA_CHECK(cudaGetLastError());      // 检查核函数启动错误
    CUDA_CHECK(cudaDeviceSynchronize()); // 等待核函数执行完成

    // 5. 设备→主机拷贝结果
    CUDA_CHECK(cudaMemcpy(out_distance, d_pairwise_dist, sizeof(float), cudaMemcpyDeviceToHost));

    // 6. 释放设备内存
    CUDA_CHECK(cudaFree(d_seq1));
    CUDA_CHECK(cudaFree(d_seq2));
    CUDA_CHECK(cudaFree(d_dtw_cost));
    CUDA_CHECK(cudaFree(d_new_dtw_cost));
    CUDA_CHECK(cudaFree(d_path_matrix));
    CUDA_CHECK(cudaFree(d_pairwise_dist));

    return 0;
}

void opendba_dtw_cleanup()
{
    cudaDeviceReset();
}