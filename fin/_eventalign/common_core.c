/* @file f5c.c
**
** f5c interface implementation
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common_core.h"
#include "common_error.h"
#include <sys/wait.h>
#include <unistd.h>

/*
todo :
Error counter for consecutive failures in the skip unreadable mode
not all the memory allocations are needed for eventalign mode
*/

/* initialise user specified options */
void init_opt(opt_t *opt)
{
    memset(opt, 0, sizeof(opt_t));
    opt->batch_size = 512;
    opt->batch_size_bases = 2 * 1000 * 1000;
    opt->pore = NULL;
    opt->num_thread = 8;
    opt->min_num_events_to_rescale = 200;
    opt->mode = 1;

#ifndef HAVE_CUDA
    opt->batch_size_bases = 5 * 1000 * 1000;
#endif
    opt->cuda_block_size = 64;
    opt->cuda_dev_id = 0;
    opt->cuda_mem_frac = 1.0f; // later set by cuda_init()

    // effective only if  CPU_GPU_PROC  is set
    opt->cuda_max_readlen = 3.0f;
    opt->cuda_avg_events_per_kmer = 2.0f; // only if CUDA_DYNAMIC_MALLOC is unset
    opt->cuda_max_avg_events_per_kmer = 5.0f;
}

/* initialise the core data structure */
core_t *init_core(opt_t opt)
{

    core_t *core = (core_t *)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    // model
    core->model = (model_t *)malloc(sizeof(model_t) * MAX_NUM_KMER);
    MALLOC_CHK(core->model);

    // load the model from files
    uint32_t kmer_size = 0;
    if (opt.mode == 1)
    {
        INFO("%s", "builtin RNA002 nucleotide model loaded");
        kmer_size = set_model(core->model, MODEL_ID_RNA002_NUCLEOTIDE);
    }
    else if (opt.mode == 2)
    {
        INFO("%s", "builtin RNA004 nucleotide model loaded");
        kmer_size = set_model(core->model, MODEL_ID_RNA004_NUCLEOTIDE);
    }
    else
    {
        ERROR("Unsupported mode %d for eventalign", opt.mode);
        exit(EXIT_FAILURE);
    }
    core->kmer_size = kmer_size;

    core->opt = opt;
    // cuda stuff
#ifdef HAVE_CUDA
    init_cuda(core);
#endif
    core->sum_bases = 0;
    // eventalign related
    core->mode = opt.mode;
    return core;
}

/* free the core data structure */
void free_core(core_t *core, opt_t opt)
{
    free(core->model);
#ifdef HAVE_CUDA
    free_cuda(core);
#endif
    free(core);
}

/* initialise a data batch */
db_t *init_db(core_t *core, int32_t read_number, int32_t ref_number)
{
    int32_t i;
    int32_t total_entries = read_number * ref_number;

    db_t *db = (db_t *)(malloc(sizeof(db_t)));
    if (db == NULL) {
        malloc_chk(NULL, __func__, __FILE__, __LINE__ - 2);
        return NULL;
    }

    db->read_id = (char **)malloc(sizeof(char *) * read_number);
    if (db->read_id == NULL) goto cleanup_db;
    db->ref_sequence = (char **)malloc(sizeof(char *) * ref_number);
    if (db->ref_sequence == NULL) goto cleanup_read_id;
    db->ref_name = (char **)malloc(sizeof(char *) * ref_number);
    if (db->ref_name == NULL) goto cleanup_ref_sequence;
    db->ref_length = (int32_t *)malloc(sizeof(int32_t) * ref_number);
    if (db->ref_length == NULL) goto cleanup_ref_name;

    db->sig = (signal_t **)malloc(sizeof(signal_t *) * read_number);
    if (db->sig == NULL) goto cleanup_ref_length;
    db->et = (event_table *)malloc(sizeof(event_table) * total_entries);
    if (db->et == NULL) goto cleanup_sig;

    db->scalings =
        (scalings_t *)malloc(sizeof(scalings_t) * total_entries);
    if (db->scalings == NULL) goto cleanup_et;

    db->event_align_pairs =
        (AlignedPair **)malloc(sizeof(AlignedPair *) * total_entries);
    if (db->event_align_pairs == NULL) goto cleanup_scalings;
    db->n_event_align_pairs =
        (int32_t *)malloc(sizeof(int32_t) * total_entries);
    if (db->n_event_align_pairs == NULL) goto cleanup_event_align_pairs;

    db->event_alignment = (event_alignment_t **)malloc(
        sizeof(event_alignment_t *) * total_entries);
    if (db->event_alignment == NULL) goto cleanup_n_event_align_pairs;
    db->n_event_alignment =
        (int32_t *)malloc(sizeof(int32_t) * total_entries);
    if (db->n_event_alignment == NULL) goto cleanup_event_alignment;

    db->events_per_base =
        (double *)malloc(sizeof(double) * total_entries);
    if (db->events_per_base == NULL) goto cleanup_n_event_alignment;

    db->base_to_event_map =
        (index_pair_t **)malloc(sizeof(index_pair_t *) * total_entries);
    if (db->base_to_event_map == NULL) goto cleanup_events_per_base;

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * total_entries);
    if (db->read_stat_flag == NULL) goto cleanup_base_to_event_map;

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * total_entries);
    if (db->site_score_map == NULL) goto cleanup_read_stat_flag;

    for (i = 0; i < total_entries; ++i)
    {
        db->site_score_map[i] = new std::map<int, ScoredSite>;
        if (db->site_score_map[i] == NULL) {
            // Clean up already allocated maps
            for (int32_t j = 0; j < i; j++) {
                delete db->site_score_map[j];
            }
            goto cleanup_site_score_map_ptr;
        }
    }

    db->total_reads = 0;

    // eventalign related

    db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * total_entries);
    if (db->event_alignment_result == NULL) goto cleanup_site_score_maps;

    db->event_alignment_result_str = (char **)malloc(sizeof(char *) * total_entries);
    if (db->event_alignment_result_str == NULL) goto cleanup_event_alignment_result;

    for (i = 0; i < total_entries; ++i)
    {
        db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        if (db->event_alignment_result[i] == NULL) {
            // Clean up already allocated vectors
            for (int32_t j = 0; j < i; j++) {
                delete db->event_alignment_result[j];
            }
            goto cleanup_event_alignment_result_str;
        }
        (db->eventalign_summary[i]).num_events = 0; // done here in the same loop for efficiency
        db->event_alignment_result_str[i] = NULL;
    }

    return db;

    // Error cleanup paths (in reverse order of allocation)
cleanup_event_alignment_result_str:
    free(db->event_alignment_result_str);
cleanup_event_alignment_result:
    free(db->event_alignment_result);
cleanup_site_score_maps:
    for (i = 0; i < total_entries; i++) {
        delete db->site_score_map[i];
    }
cleanup_site_score_map_ptr:
    free(db->site_score_map);
cleanup_read_stat_flag:
    free(db->read_stat_flag);
cleanup_base_to_event_map:
    free(db->base_to_event_map);
cleanup_events_per_base:
    free(db->events_per_base);
cleanup_n_event_alignment:
    free(db->n_event_alignment);
cleanup_event_alignment:
    free(db->event_alignment);
cleanup_n_event_align_pairs:
    free(db->n_event_align_pairs);
cleanup_event_align_pairs:
    free(db->event_align_pairs);
cleanup_scalings:
    free(db->scalings);
cleanup_et:
    free(db->et);
cleanup_sig:
    free(db->sig);
cleanup_ref_length:
    free(db->ref_length);
cleanup_ref_name:
    free(db->ref_name);
cleanup_ref_sequence:
    free(db->ref_sequence);
cleanup_read_id:
    free(db->read_id);
cleanup_db:
    free(db);
    MALLOC_CHK(NULL);
    return NULL;
}

/* load a data batch from disk */
ret_status_t load_db(core_t *core, db_t *db, int32_t read_number, int32_t ref_number,
                     char **read_id, signal_t **sig, char **ref_sequence, char **ref_name, int32_t *ref_length)
{
    double load_start = realtime();
    ret_status_t status = {0, 0};
    db->n_read_size = read_number;
    db->ref_num = ref_number;
    // set the read ids
    for (int32_t i = 0; i < read_number; ++i)
    {
        db->read_id[i] = read_id[i];
        db->sig[i] = sig[i];
    }
    // set the reference sequences
    for (int32_t i = 0; i < ref_number; ++i)
    {
        db->ref_sequence[i] = ref_sequence[i];
        db->ref_name[i] = ref_name[i];
        db->ref_length[i] = ref_length[i];
    }
    double load_end = realtime();
    INFO("Loaded batch of %d reads and %d references in %.2f sec", read_number, ref_number,
         load_end - load_start);
    status.num_reads = read_number;
    return status;
}

void event_single(core_t *core, db_t *db, int32_t i, int32_t ref_num)
{

    if (db->sig[i]->nsample > 0)
    {

        float *rawptr = db->sig[i]->rawptr;
        float range = db->sig[i]->range;
        float digitisation = db->sig[i]->digitisation;
        float offset = db->sig[i]->offset;
        int32_t nsample = db->sig[i]->nsample;

        // convert to pA
        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++)
        {
            rawptr[j] = (rawptr[j] + offset) * raw_unit;
        }

        int8_t rna = 1;
        db->et[i] = getevents(db->sig[i]->nsample, rawptr, rna);

        // get the scalings
        db->scalings[i] = estimate_scalings_using_mom(
            db->read[i], db->read_len[i], core->model, core->kmer_size, db->et[i]);

        // If sequencing RNA, reverse the events to be 3'->5'
        if (rna)
        {
            event_t *events = db->et[i].event;
            size_t n_events = db->et[i].n;
            for (size_t i = 0; i < n_events / 2; ++i)
            {
                event_t tmp_event = events[i];
                events[i] = events[n_events - 1 - i];
                events[n_events - 1 - i] = tmp_event;
            }
        }

        // allocate memory for the next alignment step
        for (int32_t j = 0; j < ref_num; ++j)
        {
            db->event_align_pairs[i * ref_num + j] = (AlignedPair *)malloc(
                sizeof(AlignedPair) * (db->et[i].n + db->ref_length[j]));
            MALLOC_CHK(db->event_align_pairs[i * ref_num + j]);
        }
    }
    else
    {
        db->et[i].n = 0;
        db->et[i].event = NULL;
        for (int32_t j = 0; j < ref_num; ++j)
        {
            db->event_align_pairs[i * ref_num + j] = NULL;
        }
    }
}

void scaling_single(core_t *core, db_t *db, int32_t i)
{

    db->event_alignment[i] = NULL;
    db->n_event_alignment[i] = 0;
    db->events_per_base[i] = 0; // todo : is double needed? not just int8?

    int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;

    if (db->n_event_align_pairs[i] > 0)
    {

        db->base_to_event_map[i] = (index_pair_t *)(malloc(sizeof(index_pair_t) * n_kmers));
        MALLOC_CHK(db->base_to_event_map[i]);

        // prepare data structures for the final calibration

        db->event_alignment[i] = (event_alignment_t *)malloc(
            sizeof(event_alignment_t) * db->n_event_align_pairs[i]);
        MALLOC_CHK(db->event_alignment[i]);

        // for (int j = 0; j < n_event_align_pairs; ++j) {
        //     fprintf(stderr, "%d-%d\n",event_align_pairs[j].ref_pos,event_align_pairs[j].read_pos);
        // }

        // todo : verify if this n is needed is needed
        db->n_event_alignment[i] = postalign(
            db->event_alignment[i], db->base_to_event_map[i], &db->events_per_base[i], db->read[i],
            n_kmers, db->event_align_pairs[i], db->n_event_align_pairs[i], core->kmer_size);

        // fprintf(stderr,"n_event_alignment %d\n",n_events);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.

        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(
            core->model, core->kmer_size, db->et[i], &db->scalings[i],
            db->event_alignment[i], db->n_event_alignment[i], 1, core->opt.min_num_events_to_rescale);

        // QC calibration
        if (!calibrated || db->scalings[i].var > MIN_CALIBRATION_VAR)
        {
            //     events[strand_idx].clear();
            free(db->event_alignment[i]);
            // free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);
    }
    else
    {
        db->base_to_event_map[i] = NULL;
        // Could not align, fail this read
        // this->events[strand_idx].clear();
        // this->events_per_base[strand_idx] = 0.0f;
        // free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_ALIGNMENT;
        return;
    }

    // Filter poor quality reads that have too many "stays"

    if (db->events_per_base[i] > 5.0)
    {
        //     events[0].clear();
        //     events[1].clear();
        // free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_QUALITY_CHK;
        return;
    }
}

/* align a single read specified by index i (perform ABEA for a single read) */
// note that this is used in f5c.cu and thus modifications must be done with care
void align_single(core_t *core, db_t *db, int32_t i)
{

    if (db->sig[i]->nsample > 0)
    { // if a good read
        if (db->sig[i]->nsample && (db->et[i].n) / (float)(db->read_len[i]) < AVG_EVENTS_PER_KMER_MAX)
        {
            db->n_event_align_pairs[i] = align(
                db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
                core->model, core->kmer_size, db->scalings[i], db->sig[i]->sample_rate);
            // fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
        }
        else
        { // todo : too many avg events per base - oversegmented
            db->n_event_align_pairs[i] = 0;
            if (core->opt.verbosity > 0)
            {
                STDERR("Skipping over-segmented read %s with %f events per base", bam_get_qname(db->bam_rec[i]), (db->et[i].n) / (float)(db->read_len[i]));
            }
        }
    }
    else
    { // if a bad read (corrupted/missing slow5 record)
        db->n_event_align_pairs[i] = 0;
    }
}
void pthread_db(core_t *core, db_t *db, void (*func)(core_t *, db_t *, int))
{

    if (core->opt.num_thread == 1)
    {
        int i;
        for (i = 0; i < db->n_bam_rec; i++)
        {
            func(core, db, i);
        }
    }
    else
    {
        // create threads
        pthread_t tids[core->opt.num_thread];
        pthread_arg_t pt_args[core->opt.num_thread];
        int32_t t, ret;
        int32_t i = 0;
        int32_t num_thread = core->opt.num_thread;
        int32_t step = (db->n_bam_rec + num_thread - 1) / num_thread;
        // todo : check for higher num of threads than the data
        // current works but many threads are created despite

        // set the data structures
        for (t = 0; t < num_thread; t++)
        {
            pt_args[t].core = core;
            pt_args[t].db = db;
            pt_args[t].starti = i;
            i += step;
            if (i > db->n_bam_rec)
            {
                pt_args[t].endi = db->n_bam_rec;
            }
            else
            {
                pt_args[t].endi = i;
            }
            pt_args[t].func = func;
#ifdef WORK_STEAL
            pt_args[t].all_pthread_args = (void *)pt_args;
#endif
            // fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);
        }

        // create threads
        for (t = 0; t < core->opt.num_thread; t++)
        {
            ret = pthread_create(&tids[t], NULL, pthread_single,
                                 (void *)(&pt_args[t]));
            NEG_CHK(ret);
        }

        // pthread joining
        for (t = 0; t < core->opt.num_thread; t++)
        {
            int ret = pthread_join(tids[t], NULL);
            NEG_CHK(ret);
        }
    }
}
/* align a data batch (perform ABEA for a data batch) */
void align_db(core_t *core, db_t *db)
{
#ifdef HAVE_CUDA
    align_cuda(core, db);
#else
    pthread_db(core, db, align_single);
#endif
}

void eventalign_single(core_t *core, db_t *db, int32_t i)
{
    realign_read(db->event_alignment_result[i], &(db->eventalign_summary[i]), core->event_summary_fp, db->fasta_cache[i], core->m_hdr,
                 db->bam_rec[i], db->read_len[i],
                 i,
                 core->clip_start,
                 core->clip_end,
                 &(db->et[i]), core->model, core->kmer_size, db->base_to_event_map[i], db->scalings[i], db->events_per_base[i], db->sig[i]->sample_rate);

    char *qname = bam_get_qname(db->bam_rec[i]);
    char *contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];
    std::vector<event_alignment_t> *event_alignment_result = db->event_alignment_result[i];
    int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0;
    int8_t scale_events = (core->opt.flag & F5C_SCALE_EVENTS) ? 1 : 0;
    int8_t collapse_events = (core->opt.flag & F5C_COLLAPSE_EVENTS) ? 1 : 0;
    int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0;
    int8_t write_signal_index = (core->opt.flag & F5C_PRINT_SIGNAL_INDEX) ? 1 : 0;
    int8_t sam_output = (core->opt.flag & F5C_SAM) ? 1 : 0;
    int8_t paf_output = (core->opt.flag & F5C_PAF) ? 1 : 0;
    int8_t m6anet_output = (core->opt.flag & F5C_M6ANET) ? 1 : 0;
    int8_t rna = (core->opt.flag & F5C_RNA) ? 1 : 0;

    if (paf_output)
    {
        int64_t ref_len = core->m_hdr->target_len[db->bam_rec[i]->core.tid];
        db->event_alignment_result_str[i] = emit_event_alignment_paf(&(db->et[i]), db->sig[i]->nsample, ref_len, core->kmer_size, db->scalings[i], *event_alignment_result, db->bam_rec[i], qname, contig, rna);
    }
    else if (sam_output)
    {
        int8_t sam_out_version = core->opt.sam_out_version;
        int64_t ref_len = core->m_hdr->target_len[db->bam_rec[i]->core.tid];
        db->event_alignment_result_str[i] = emit_event_alignment_sam(qname, core->m_hdr, db->bam_rec[i], *event_alignment_result, sam_out_version, &(db->et[i]), db->sig[i]->nsample, ref_len, rna, db->scalings[i]);
    }
    else if (m6anet_output)
    {
        db->event_alignment_result_str[i] = emit_event_alignment_tsv_m6anet(0, &(db->et[i]), core->model, core->kmer_size, db->scalings[i], *event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index, collapse_events,
                                                                            db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);
    }
    else
    {
        db->event_alignment_result_str[i] = emit_event_alignment_tsv(0, &(db->et[i]), core->model, core->kmer_size, db->scalings[i], *event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index, collapse_events,
                                                                     db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);
    }
}

void process_single(core_t *core, db_t *db, int32_t i)
{
    event_single(core, db, i);
    align_single(core, db, i);
    scaling_single(core, db, i);
}

/* completely process a data batch
   (all steps: event detection, adaptive banded event alignment, ...., HMM) */
void process_db(core_t *core, db_t *db)
{
    double process_start = realtime();

    if ((core->opt.flag & F5C_SEC_PROF) || (!(core->opt.flag & F5C_DISABLE_CUDA)))
    {

        double realtime0 = core->realtime0;

        double event_start = realtime();
        pthread_db(core, db, event_single);
        double event_end = realtime();
        core->event_time += (event_end - event_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double align_start = realtime();
        align_db(core, db);
        double align_end = realtime();
        core->align_time += (align_end - align_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double est_scale_start = realtime();
        pthread_db(core, db, scaling_single);
        double est_scale_end = realtime();
        core->est_scale_time += (est_scale_end - est_scale_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double meth_start = realtime();
        pthread_db(core, db, meth_single);
        double meth_end = realtime();
        core->meth_time += (meth_end - meth_start);

        fprintf(stderr, "[%s::%.3f*%.2f] HMM done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }
    else
    {
        if (core->opt.num_thread == 1)
        {
            int32_t i = 0;
            for (i = 0; i < db->n_bam_rec; i++)
            {
                process_single(core, db, i);
            }
        }
        else
        {
            pthread_db(core, db, process_single);
        }
    }

    double process_end = realtime();
    core->process_db_time += (process_end - process_start);

    return;
}

/* write the output for a processed data batch */
void output_db(core_t *core, db_t *db)
{

    double output_start = realtime();

    if (core->opt.flag & F5C_PRINT_EVENTS)
    {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++)
        {
            printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
                   bam_get_qname(db->bam_rec[i]), (int)db->et[i].n,
                   (int)db->et[i].start, (int)db->et[i].end);
            uint32_t j = 0;
            for (j = 0; j < db->et[i].n; j++)
            {
                printf("{%d,%f,%f,%f}\t", (int)db->et[i].event[j].start,
                       db->et[i].event[j].length, db->et[i].event[j].mean,
                       db->et[i].event[j].stdv);
            }
            printf("\n");
        }
    }
    if (core->opt.flag & F5C_PRINT_BANDED_ALN)
    {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++)
        {
            if ((db->read_stat_flag[i]) & FAILED_ALIGNMENT)
            {
                continue;
            }
            printf(">%s\tN_ALGN_PAIR:%d\t{ref_pos,read_pos}\n",
                   bam_get_qname(db->bam_rec[i]),
                   (int)db->n_event_align_pairs[i]);
            AlignedPair *event_align_pairs = db->event_align_pairs[i];
            int32_t j = 0;
            for (j = 0; j < db->n_event_align_pairs[i]; j++)
            {
                printf("{%d,%d}\t", event_align_pairs[j].ref_pos,
                       event_align_pairs[j].read_pos);
            }
            printf("\n");
        }
    }

    if (core->opt.flag & F5C_PRINT_SCALING)
    {
        int32_t i = 0;
        printf("read\tshift\tscale\tvar\n");

        for (i = 0; i < db->n_bam_rec; i++)
        {
            if ((db->read_stat_flag[i]) & (FAILED_ALIGNMENT | FAILED_CALIBRATION))
            {
                continue;
            }
            printf("%s\t%.2lf\t%.2lf\t%.2lf\n", bam_get_qname(db->bam_rec[i]),
                   db->scalings[i].shift, db->scalings[i].scale,
                   db->scalings[i].var);
        }
    }

    core->sum_bases += db->sum_bases;
    core->total_reads += db->total_reads;
    core->bad_fast5_file += db->bad_fast5_file;
    core->ultra_long_skipped += db->ultra_long_skipped;
    core->skip_mapq_reads += db->skip_mapq_reads;
    core->skip_sec_reads += db->skip_sec_reads;
    core->unmapped_reads += db->unmapped_reads;

    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; i++)
    {
        if (!db->read_stat_flag[i])
        {
            char *qname = bam_get_qname(db->bam_rec[i]);
            char *contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];

            if (core->mode == 0)
            {
                std::map<int, ScoredSite> *site_score_map = db->site_score_map[i];
                // write all sites for this read
                for (auto iter = site_score_map->begin(); iter != site_score_map->end(); ++iter)
                {

                    const ScoredSite &ss = iter->second;
                    double sum_ll_m = ss.ll_methylated[0];   //+ ss.ll_methylated[1];
                    double sum_ll_u = ss.ll_unmethylated[0]; //+ ss.ll_unmethylated[1];
                    double diff = sum_ll_m - sum_ll_u;

                    // output only if inside the window boundaries
                    if (!((core->clip_start != -1 && ss.start_position < core->clip_start) ||
                          (core->clip_end != -1 && ss.end_position >= core->clip_end)))
                    {
                        if (core->opt.meth_out_version == 1)
                        {
                            printf("%s\t%d\t%d\t", contig, ss.start_position, ss.end_position);
                        }
                        else if (core->opt.meth_out_version == 2)
                        {
                            printf("%s\t%c\t%d\t%d\t", contig, bam_is_rev(db->bam_rec[i]) ? '-' : '+', ss.start_position, ss.end_position);
                        }
                        printf("%s\t%.2lf\t", qname, diff);
                        printf("%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
                        printf("%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());
                    }
                }
            }

            else if (core->mode == 1)
            {
                FILE *summary_fp = core->event_summary_fp;
                EventalignSummary summary = db->eventalign_summary[i];
                scalings_t scalings = db->scalings[i];
                if (summary_fp != NULL && summary.num_events > 0)
                {
                    size_t strand_idx = 0;
                    std::string fast5_path_str = core->readbb->get_signal_path(qname);
                    const char *path = (core->opt.flag & F5C_RD_SLOW5) ? "slow5" : fast5_path_str.c_str();
                    fprintf(summary_fp, "%ld\t%s\t", (long)(db->read_idx[i]), qname);
                    fprintf(summary_fp, "%s\t%s\t%s\t", path, (core->opt.flag & F5C_RNA) ? "rna" : "dna", strand_idx == 0 ? "template" : "complement");
                    fprintf(summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
                    fprintf(summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration / (db->sig[i]->sample_rate), scalings.shift, scalings.scale, 0.0, scalings.var);
                }

                char *event_alignment_result_str = db->event_alignment_result_str[i];
                fputs(event_alignment_result_str, stdout);
            }
        }
        else
        {
            if ((db->read_stat_flag[i]) & FAILED_CALIBRATION)
            {
                core->failed_calibration_reads++;
            }
            else if ((db->read_stat_flag[i]) & FAILED_ALIGNMENT)
            {
                core->failed_alignment_reads++;
            }
            else if ((db->read_stat_flag[i]) & FAILED_QUALITY_CHK)
            {
                core->qc_fail_reads++;
            }
            else
            {
                assert(0);
            }
        }
    }
    fflush(stdout);

    // core->read_index = core->read_index + db->n_bam_rec;
    double output_end = realtime();
    core->output_time += (output_end - output_start);
}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t *db)
{
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i)
    {
        bam_destroy1(db->bam_rec[i]);
        db->bam_rec[i] = bam_init1();
        free(db->fasta_cache[i]);
        free(db->read[i]);
        free(db->sig[i]->rawptr);
        free(db->sig[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
        delete db->site_score_map[i];
        db->site_score_map[i] = new std::map<int, ScoredSite>;

        if (db->event_alignment_result)
        { // eventalign related
            delete db->event_alignment_result[i];
            db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        }
        if (db->event_alignment_result_str)
        { // eventalign related
            free(db->event_alignment_result_str[i]);
            db->event_alignment_result_str[i] = NULL;
        }
    }
}

/* completely free a data batch */
void free_db(db_t *db)
{
    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i)
    {
        bam_destroy1(db->bam_rec[i]);
    }
    free(db->bam_rec);
    free(db->fasta_cache);
    free(db->read);
    free(db->read_len);
    free(db->read_idx);
    free(db->et);
    free(db->sig);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);
    for (i = 0; i < db->capacity_bam_rec; ++i)
    {
        delete db->site_score_map[i];
    }
    free(db->site_score_map);
    // eventalign related
    if (db->eventalign_summary)
    {
        free(db->eventalign_summary);
    }
    if (db->event_alignment_result)
    {
        for (i = 0; i < db->capacity_bam_rec; ++i)
        {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }
    if (db->event_alignment_result_str)
    {
        free(db->event_alignment_result_str);
    }

    free(db);
}
