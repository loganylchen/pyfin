
/* @file model.c
**
** load a pore model from file or memory
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "model.h"
#include "error.h"
#include "common.h"

uint32_t eval_num_kmer(uint32_t kmer_size, uint32_t type)
{
    uint32_t num_kmer = 0;
    num_kmer = (uint32_t)(1 << 2 * kmer_size); // num_kmer should be 4^kmer_size
    assert(num_kmer <= MAX_NUM_KMER);
    return num_kmer;
}

uint32_t set_model(model_t *model, uint32_t model_id)
{

    uint32_t kmer_size = 0;
    uint32_t num_kmer = 0;
    float *inbuilt_model = NULL;

    if (model_id == MODEL_ID_RNA002_NUCLEOTIDE)
    {
        kmer_size = 5;
        num_kmer = 1024;
        inbuilt_model = rna002_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2 * kmer_size)); // num_kmer should be 4^kmer_size
    }
    else if (model_id == MODEL_ID_RNA004_NUCLEOTIDE)
    {
        kmer_size = 9;
        num_kmer = 262144;
        inbuilt_model = rna004_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2 * kmer_size)); // num_kmer should be 4^kmer_size
    }
    else
    {
        assert(0);
    }

    uint32_t i = 0;
    for (i = 0; i < num_kmer; i++)
    {
        model[i].level_mean = inbuilt_model[i * 2 + 0];
        model[i].level_stdv = inbuilt_model[i * 2 + 1];
#ifdef CACHED_LOG
        model[i].level_log_stdv = log(model[i].level_stdv);
#endif
    }

#ifdef DEBUG_MODEL_PRINT
    i = 0;
    fprintf(stderr, "level_mean\tlevel_stdv\tsd_mean\tsd_stdv\n");
    for (i = 0; i < num_kmer; i++)
    {
        fprintf(stderr, "%f\t%f\t%f\t%f\n", model[i].level_mean,
                model[i].level_stdv, 0.0, 0.0);
    }
#endif

    return kmer_size;
}
