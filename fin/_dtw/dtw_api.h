#ifndef DTW_C_API_H
#define DTW_C_API_H

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * @brief 调用 OpenDBA 原版 CUDA DTW 计算两个序列的距离
     * @param seq1 主机端浮点序列1（float 类型）
     * @param len1 序列1长度
     * @param seq2 主机端浮点序列2（float 类型）
     * @param len2 序列2长度
     * @param use_open_start 是否启用 open start 边界
     * @param use_open_end 是否启用 open end 边界
     * @param out_distance 输出 DTW 距离（主机端，float 类型）
     * @return 0=成功，非0=错误码（1=内存分配失败，2=核函数启动失败，3=数据拷贝失败）
     */
    int opendba_dtw_cuda(
        const float *seq1, size_t len1,
        const float *seq2, size_t len2,
        int use_open_start,
        int use_open_end,
        float *out_distance);

    /**
     * @brief 释放 CUDA 资源
     */
    void opendba_dtw_cleanup();

#ifdef __cplusplus
}
#endif

#endif // DTW_C_API_H