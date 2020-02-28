#include "Compressor.h"
#include <omp.h>

#ifdef __DEBUG_LOG__
FILE *fp_log_debug = NULL;
unsigned int debug_count = 0;
#endif

#ifdef __GLOBAL_STATS__
double under_T = 0;
double over_T = 0;
#endif

/* --------------------------------------------------------------------------
 * Compression functions.
 */

char* in_buf;
uint AVG_CANT, B_CTX, B_CTX_LEN, NS_MODEL_SIZE;

uint16_t* ctx_avgs_sums;
uint16_t* ctx_avgs_err_sums;
uint32_t* ctx_err_avgs_total;

void init_global_stats(context_models* cm) {

    cm->model_seq8 = new BASE_MODEL<uint8_t>[NS_MODEL_SIZE];
    ctx_avgs_sums = new uint16_t[AVG_CANT];
    ctx_avgs_err_sums = new uint16_t[AVG_CANT];
    ctx_err_avgs_total = new uint32_t[Q_CTX];

    memset(ctx_avgs_sums, 0, AVG_CANT * sizeof(uint16_t));
    memset(ctx_avgs_err_sums, 0 << AVG_SHIFT, AVG_CANT * sizeof(uint16_t));

    for (uint s_ctx = 0; s_ctx < B_CTX; s_ctx++) {
        for (uint dif = 0; dif < DIF_CANT; dif++) {
            for (uint q_quant = 0; q_quant < Q_LOG_CANT; q_quant++) {
                uint avg_ctx = (s_ctx << TOTAL_Q_LOG) + (dif << Q_LOG) + q_quant;
                ctx_avgs_sums[avg_ctx] = q_quant << AVG_SHIFT;
            }
        }
    }

    memset(ctx_err_avgs_total, 1 << TOTAL_ERR_SHIFT, Q_CTX * sizeof(uint32_t));
}

void delete_global_stats(context_models* cm) {
    delete [] cm->model_seq8;
    delete cm;
    delete [] ctx_avgs_sums;
    delete [] ctx_avgs_err_sums;
    delete [] ctx_err_avgs_total;
}

void copy_average_stats(Compressor* c) {
    memcpy(c->ctx_avgs_sums, ctx_avgs_sums, AVG_CANT* sizeof(uint16_t));
    memcpy(c->ctx_avgs_err_sums, ctx_avgs_err_sums, AVG_CANT* sizeof(uint16_t));
    memcpy(c->ctx_err_avgs_total, ctx_err_avgs_total, Q_CTX * sizeof(uint32_t));
}

/*
 * A blocking read that refuses to return truncated reads.
 */
static ssize_t xread(int fd, char *buf, size_t count) {
    ssize_t len, tlen;
//    printf("Read %d \n", read_num++);
    tlen = 0;
    do {
        len = read(fd, buf, count);
        if (len == -1) {
            if (errno == EINTR)
                continue;
            return -1;
        }

        if (len == 0)
            return tlen;

        buf += len;
        count -= len;
        tlen += len;
    } while (count);

    return tlen;
}

bool load_data(int in_fd, Compressor ** comps, int update_load, uint &blocks_loaded, int &blk_start) {

    int sz;

    u_char comp_id = 0;
    blocks_loaded = 0;

    while (comp_id < update_load && (sz = xread(in_fd, &in_buf[blk_start], BLK_SIZE - blk_start)) > 0) {

        char *in_end = NULL;
        int nseqs = 0, remainder_length = 0, end = 0;

        Compressor* c = comps[comp_id];

        int error = 0;
        if (0 > (error = c->fq_parse_reads(in_buf, sz + blk_start,
                                    &in_end, &nseqs, &remainder_length, end))) {
            fprintf(stderr, "Failure to parse and/or compress. Error %d \n", error);
            return true;
        }

        blocks_loaded++;
        comp_id++;

        /* We maybe ended on a partial fastq entry, so start from there */
        if (remainder_length > 0) {
            memmove(in_buf, in_end, remainder_length);
            blk_start = remainder_length - 1;
        } else
            blk_start = 0;
    }

    return blocks_loaded <= 0;
}

void update_stats(context_models* cm, Compressor** comps, uc blocks_loaded){

    uint i;

    for (i = 0; i < AVG_CANT; i++) {
        uint sum_ctx_avgs_sums = 0;
        uint sum_ctx_avgs_err_sums = 0;
        for (uc c = 0; c < blocks_loaded; c++) {
            sum_ctx_avgs_sums += comps[c]->ctx_avgs_sums[i];
            sum_ctx_avgs_err_sums += comps[c]->ctx_avgs_err_sums[i];
        }
        ctx_avgs_sums[i] = round((double)sum_ctx_avgs_sums/blocks_loaded);
        ctx_avgs_err_sums[i] = round((double)sum_ctx_avgs_err_sums/blocks_loaded);
    }

    for (i = 0; i < Q_CTX; i++) {
        uint sum_ctx_err_avgs_total = 0;
        for (uc c = 0; c < blocks_loaded; c++) {
            sum_ctx_err_avgs_total += comps[c]->ctx_err_avgs_total[i];
        }
        ctx_err_avgs_total[i] = round((double)sum_ctx_err_avgs_total/blocks_loaded);
    }

    void **models = new void *[blocks_loaded];

    for (i = 0; i < NS_MODEL_SIZE; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_seq8[i]);
        }
        cm->model_seq8[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < 256; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_name_prefix[i]);
        }
        cm->model_name_prefix[i].mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_name_suffix[i]);
        }
        cm->model_name_suffix[i].mix_array(models, blocks_loaded);

        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_name_len[i]);
        }
        cm->model_name_len[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < 8192; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_name_middle[i]);
        }
        cm->model_name_middle[i].mix_array(models, blocks_loaded);
    }

    for (i = 0; i < CTX_CNT; i++) {
        for (uint c = 0; c < blocks_loaded; c++) {
            models[c] = (void*)(&comps[c]->cm->model_qual_quant[i]);
        }
        cm->model_qual_quant[i].mix_array(models, blocks_loaded);
    }

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void*)(&comps[c]->cm->quant_top);
    }
    cm->quant_top.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void*)(&comps[c]->cm->model_len1);
    }
    cm->model_len1.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void*)(&comps[c]->cm->model_len2);
    }
    cm->model_len2.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void*)(&comps[c]->cm->model_len3);
    }
    cm->model_len3.mix_array(models, blocks_loaded);

    for (uint c = 0; c < blocks_loaded; c++) {
        models[c] = (void*)(&comps[c]->cm->model_same_len);
    }
    cm->model_same_len.mix_array(models, blocks_loaded);

    delete [] models;
}

/*
 * Encode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */
/*
* Parse one block at a time. Blocks may not terminate on exact fastq
* boundaries, so we need to know where we ended processing and move
* that partial fastq entry back to the start for the next block.
*
* We write out the block size too so we can decompress block at a time.
*/
int encode(int in_fd, int out_fd, enano_params* p) {

    int res = 0;
    bool decode = false;
    double start_time = omp_get_wtime();
    double enc_time = 0, code_time = 0, load_time = 0, write_time = 0, update_time = 0;

    printf("Starting encoding in FAST MODE with %d threads... \n", p->num_threads);

    int blk_start = 0;
    uint block_num = 0;

    uint BLK_UPD_FREQ = p->blk_upd_freq;
    uint BLK_UPD_THRESH = p->blk_upd_thresh + 1;

    u_char cant_compressors = MAX(BLK_UPD_FREQ, p->num_threads);

    Compressor** comps = new Compressor*[cant_compressors];
    for (uint i = 0; i < cant_compressors; i ++) {
        comps[i] = new Compressor(p);
    }

    printf("Starting adaptative encoding for %d (%d + 1) blocks, and update every %d blocks... \n", BLK_UPD_THRESH, BLK_UPD_THRESH - 1, BLK_UPD_FREQ);

    bool finished = false;

    context_models* cm = new context_models;
    init_global_stats(cm);

    //Update stats
    uint blocks_loaded;
    uint update_blocks = 0;

    uint num_batches = 1 + BLK_UPD_THRESH/BLK_UPD_FREQ;
    uint* update_batches = new uint[1 + BLK_UPD_THRESH/BLK_UPD_FREQ];
    update_batches[0] = 1;
    for (uint i = 1; i < num_batches; i++)
        update_batches[i] = BLK_UPD_FREQ;

    uint batch = 0;
    uint update_load = update_batches[batch];//MIN(BLK_UPD_FREQ, BLK_UPD_THRESH);

    while (update_blocks < BLK_UPD_THRESH && !(finished = load_data(in_fd, comps, update_load,blocks_loaded,blk_start))) {
        #pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            comps[i]->soft_reset();
            copy_average_stats(comps[i]);
            comps[i]->copy_stats(cm);
            comps[i]->fq_compress();
        }

        for (uint i = 0; i < blocks_loaded; i++) {
            if (!comps[i]->output_block(out_fd)) {
                fprintf(stderr, "Abort: truncated write.\n");
                finished = true;
                res = -1;
                break;
            }
        }
        update_stats(cm, comps, blocks_loaded);

        update_blocks += blocks_loaded;
        block_num += blocks_loaded;
        batch += 1;
        update_load = MIN((BLK_UPD_THRESH - update_blocks),  update_batches[batch]);
    }

    delete [] update_batches;

    printf("Starting parallelized fast encoding...\n");
    //Update context models accumulated probabilities.
    comps[0]->update_AccFreqs(cm, decode);
    //No updates from now on
    for (uint i = 0; i < cant_compressors; i++) {
        comps[i]->updateModel = false;
        delete [] comps[i]->cm->model_seq8;
        delete comps[i]->cm;
        comps[i]->cm = cm;
    }

    //Finished updating the models
    update_time += omp_get_wtime() - start_time;

    //Parallelized compression with fixed stats
    double clock = 0;
    while (!finished) {
        uint blocks_loaded = 0;

        //Load and parse the data for each compressor
        clock = omp_get_wtime();
        finished = load_data(in_fd, comps, cant_compressors, blocks_loaded, blk_start);
        load_time += omp_get_wtime() - clock;
        //Encode the loaded data in parallel in each of the threads
        clock = omp_get_wtime();
        #pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            comps[i]->soft_reset();
            copy_average_stats(comps[i]);
            comps[i]->fq_compress();
        }
        code_time += omp_get_wtime() - clock;
        //Write the data to the file in order
        clock = omp_get_wtime();
        for (uint i = 0; i < blocks_loaded; i++) {
            if (!comps[i]->output_block(out_fd)) {
                fprintf(stderr, "Abort: truncated write.\n");
                finished = true;
                res = -1;
                break;
            }
        }
        write_time += omp_get_wtime() - clock;
        block_num += blocks_loaded;
    }

    printf("Total encoded blocks: %d \n", block_num);

    long long name_in = 0, name_out = 0, base_in = 0, base_out = 0, qual_in = 0, qual_out = 0;
    for (uint i = 0; i < cant_compressors; i++) {
        name_in += comps[i]->name_in;
        name_out += comps[i]->name_out;
        base_in += comps[i]->base_in;
        base_out += comps[i]->base_out;
        qual_in += comps[i]->qual_in;
        qual_out += comps[i]->qual_out;
    }

    for (uint i = 0; i < cant_compressors; i ++)
        delete comps[i];

    delete [] comps;

    delete_global_stats(cm);

    enc_time = omp_get_wtime() - start_time;

#ifdef __TIMING__
    fprintf(stdout, "Update models: %.2f s\n",
            (double)update_time);
    fprintf(stdout, "Load time: %.2f s\n",
            (double)load_time);
    fprintf(stdout, "Coding time: %.2lf s\n",
            (double)code_time);
    fprintf(stdout, "Write time: %.2f s\n",
            (double)write_time);
#endif
    
    fprintf(stdout, "Stream <original size in bytes> -> <compressed size in bytes> (<compression ratio>)\n");

    fprintf(stdout, "IDs   %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            name_in, name_out, (double) name_out / name_in);
    fprintf(stdout, "Bases %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            base_in, base_out, (double) base_out / base_in);
    fprintf(stdout, "Quals %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            qual_in, qual_out, (double) qual_out / qual_in);
    fprintf(stdout, "Total %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            name_in + base_in + qual_in, name_out + base_out + qual_out, (double) (name_out + base_out + qual_out) / (name_in + base_in + qual_in));
    fprintf(stdout, "Total compression time: %.2f s\n",
            (double)enc_time);

#ifdef __GLOBAL_STATS__
    fprintf(stdout, "Under T %0.1f\% Over T %0.1f\% \n",
            under_T * 100 / (over_T + under_T), over_T * 100 / (over_T + under_T));
#endif

    return res;
}


int encode_st(int in_fd, int out_fd, enano_params* p) {

    int res = 0;
    double start_time = omp_get_wtime();
    double enc_time = 0, code_time = 0, load_time = 0, write_time = 0, update_time = 0;

    printf("Starting encoding in Max Compression mode... \n");

    int blk_start = 0;
    uint block_num = 0;

    uint cant_compressors = 1;

    Compressor** comps = new Compressor*[cant_compressors];
    comps[0] = new Compressor(p);

    //Update stats
    uint blocks_loaded;
    while (!load_data(in_fd, comps, 1,blocks_loaded,blk_start)) {
        comps[0]->fq_compress();
        if (!comps[0]->output_block(out_fd)) {
            fprintf(stderr, "Abort: truncated write.\n");
            res = -1;
            break;
        }
        block_num += blocks_loaded;
    }

    printf("Total encoded blocks: %d \n", block_num);

    long long name_in = 0, name_out = 0, base_in = 0, base_out = 0, qual_in = 0, qual_out = 0;
    name_in += comps[0]->name_in;
    name_out += comps[0]->name_out;
    base_in += comps[0]->base_in;
    base_out += comps[0]->base_out;
    qual_in += comps[0]->qual_in;
    qual_out += comps[0]->qual_out;


    delete comps[0];
    delete [] comps;

    enc_time = omp_get_wtime() - start_time;

#ifdef __TIMING__
    fprintf(stdout, "Update models: %.2f s\n",
            (double)update_time);
    fprintf(stdout, "Load time: %.2f s\n",
            (double)load_time);
    fprintf(stdout, "Coding time: %.2lf s\n",
            (double)code_time);
    fprintf(stdout, "Write time: %.2f s\n",
            (double)write_time);
#endif
    
    fprintf(stdout, "Stream <original size in bytes> -> <compressed size in bytes> (<compression ratio>)\n");

    fprintf(stdout, "IDs   %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            name_in, name_out, (double) name_out / name_in);
    fprintf(stdout, "Bases %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            base_in, base_out, (double) base_out / base_in);
    fprintf(stdout, "Quals %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            qual_in, qual_out, (double) qual_out / qual_in);
    fprintf(stdout, "Total %" PRIu64 " -> %" PRIu64 " (%0.3f)\n",
            name_in + base_in + qual_in, name_out + base_out + qual_out, (double) (name_out + base_out + qual_out) / (name_in + base_in + qual_in));
    fprintf(stdout, "Total compression time: %.2f s\n",
            (double)enc_time);

#ifdef __GLOBAL_STATS__
    fprintf(stdout, "Under T %0.1f\% Over T %0.1f\% \n",
            under_T * 100 / (over_T + under_T), over_T * 100 / (over_T + under_T));
#endif

    return res;
}

bool load_data_decode(int in_fd, Compressor ** comps, int update_load, uint &blocks_loaded) {

    bool res = 0;
    unsigned char len_buf[4];
    int sz;

    blocks_loaded = 0;

    u_char comp_id = 0;

    while ( comp_id < update_load && (4 == (sz = read(in_fd, len_buf, 4))) ) {
        int32_t comp_len =
                (len_buf[0] << 0) +
                (len_buf[1] << 8) +
                (len_buf[2] << 16) +
                (len_buf[3] << 24);

        int rem_len = comp_len, in_off = 0;

        do {
            errno = 0;
            int tmp_len = read(in_fd, comps[comp_id]->decode_buf + in_off, rem_len);
            if (errno == EINTR && tmp_len == -1)
                continue;

            if (tmp_len == -1) {
                fprintf(stderr, "Abort: read failed, %d.\n", errno);
                perror("foo");
                res = -1;
                goto error;
            }
            if (tmp_len == 0) {
                fprintf(stderr, "Abort: truncated read, %d.\n", errno);
                res = -1;
                goto error;
            }
            rem_len -= tmp_len;
            in_off += tmp_len;
        } while (rem_len);

        blocks_loaded++;
        comp_id++;
    }
    error:
    return blocks_loaded <= 0 || res == -1;
}

/*
 * Decode an entire stream
 *
 * Returns 0 on success
 *        -1 on failure
 */

/*
 * Parse one block at a time. Blocks may not terminate on exact fastq
 * boundaries, so we need to know where we ended processing and move
 * that partial fastq entry back to the start for the next block.
 *
 * We write out the block size too so we can decompress block at a time.
 */
int decode(int in_fd, int out_fd, enano_params* p) {
    int res = 0;
    bool decode = true;
    double start_time = omp_get_wtime();
    double dec_time = 0, decode_time = 0, load_time = 0, write_time = 0, update_time = 0;
    double clock = 0;

    printf("Starting decoding with %d threads... \n", p->num_threads);

    uint block_num = 0;

    uint BLK_UPD_FREQ = p->blk_upd_freq;
    uint BLK_UPD_THRESH = p->blk_upd_thresh + 1;

    u_char cant_compressors = MAX(BLK_UPD_FREQ, p->num_threads);

    Compressor** comps = new Compressor*[cant_compressors];
    for (uint i = 0; i < cant_compressors; i ++) {
        comps[i] = new Compressor(p);
        comps[i]->decode_buf = new char[BLK_SIZE];
    }

    printf("Starting decoding with context model update... \n");

    bool finished = false;

    context_models* cm = new context_models;
    init_global_stats(cm);

    //Update stats
    uint blocks_loaded;
    uint update_blocks = 0;
    uint num_batches = 1 + BLK_UPD_THRESH/BLK_UPD_FREQ;
    uint* update_batches = new uint[1 + BLK_UPD_THRESH/BLK_UPD_FREQ];
    update_batches[0] = 1;
    for (uint i = 1; i < num_batches; i++)
        update_batches[i] = BLK_UPD_FREQ;

    uint batch = 0;
    uint update_load = update_batches[batch];

    while (update_blocks < BLK_UPD_THRESH && !(finished = load_data_decode(in_fd, comps, update_load, blocks_loaded))) {

#pragma omp parallel for
        for (uint i = 0; i < blocks_loaded; i++) {
            comps[i]->soft_reset();
            copy_average_stats(comps[i]);
            comps[i]->copy_stats(cm);
            comps[i]->fq_decompress();
        }
        //Write output
        for (uint i = 0; i < blocks_loaded; i++) {
            if (comps[i]->uncomp_len != write(out_fd, comps[i]->out_buf, comps[i]->uncomp_len)) {
                fprintf(stderr, "Abort: truncated write.\n");
                res = -1;
                goto finishdecode;
            }
        }
        update_stats(cm, comps, blocks_loaded);

        update_blocks += blocks_loaded;
        block_num += blocks_loaded;
        batch += 1;
        update_load = MIN((BLK_UPD_THRESH - update_blocks),  update_batches[batch]);
    }

    delete [] update_batches;

    if (!finished) {

        printf("Starting parallelized fast decoding... \n");

        //Update context models accumulated probabilities.
        comps[0]->update_AccFreqs(cm, decode);
        //No updates from now on
        for (uint i = 0; i < cant_compressors; i++) {
            comps[i]->updateModel = false;
            delete comps[i]->cm;
            comps[i]->cm = cm;
        }

        //Finished updating the models
        update_time += omp_get_wtime() - start_time;

        //Parallelized decompression with fixed stats
        while (!finished) {

            uint blocks_loaded = 0;

            clock = omp_get_wtime();
            finished = load_data_decode(in_fd, comps, cant_compressors, blocks_loaded);
            load_time += omp_get_wtime() - clock;

            clock = omp_get_wtime();
            #pragma omp parallel for
            for (uint i = 0; i < blocks_loaded; i ++) {
                comps[i]->soft_reset();
                copy_average_stats(comps[i]);
                comps[i]->fq_decompress();
            }
            decode_time += omp_get_wtime() - clock;

            clock = omp_get_wtime();
            for (uint i = 0; i < blocks_loaded; i ++) {
                if (comps[i]->uncomp_len != write(out_fd, comps[i]->out_buf, comps[i]->uncomp_len)) {
                    fprintf(stderr, "Abort: truncated write.\n");
                    res = -1;
                    goto finishdecode;
                }
            }
            write_time += omp_get_wtime() - clock;

            block_num += blocks_loaded;
        }
    }
    //We use this goto flag to break the double loop
    finishdecode:

    dec_time = omp_get_wtime() - start_time;

    printf("Total decoded blocks: %d\n", block_num, BLK_UPD_THRESH);

#ifdef __TIMING__
    fprintf(stdout, "Update models: %.2f s\n",
            (double)update_time);
    fprintf(stdout, "Load time: %.2f s\n",
            (double)load_time);
    fprintf(stdout, "Coding time: %.2lf s\n",
            (double)decode_time);
    fprintf(stdout, "Write time: %.2f s\n",
            (double)write_time);
#endif
    fprintf(stdout, "Total decompression time: %.2f s\n",
            (double)dec_time);

    for (uint i = 0; i < cant_compressors; i ++) {
        delete comps[i]->decode_buf;
        if (comps[i]->cm != cm)
            delete comps[i]->cm;
        delete comps[i];
    }
    delete [] comps;

    delete_global_stats(cm);

    return res;
}

int decode_st (int in_fd, int out_fd, enano_params* p) {
    int res = 0;
    double start_time = omp_get_wtime();
    double dec_time = 0, decode_time = 0, load_time = 0, write_time = 0, update_time = 0;

    printf("Starting decoding Max Compression mode... \n");

    uint block_num = 0;

    u_char cant_compressors = 1;

    Compressor** comps = new Compressor*[cant_compressors];
    comps[0] = new Compressor(p);
    comps[0]->decode_buf = new char[BLK_SIZE];

    //Update stats
    uint blocks_loaded;
    while (!load_data_decode(in_fd, comps, 1, blocks_loaded)) {
        comps[0]->fq_decompress();
        //Write output
        if (comps[0]->uncomp_len != write(out_fd, comps[0]->out_buf, comps[0]->uncomp_len)) {
            fprintf(stderr, "Abort: truncated write.\n");
            res = -1;
            break;
        }
        block_num += blocks_loaded;
    }

    delete comps[0]->decode_buf;
    delete comps[0];
    delete [] comps;

    dec_time = omp_get_wtime() - start_time;

    printf("Total decoded blocks: %d\n", block_num);

#ifdef __TIMING__
    fprintf(stdout, "Update models: %.2f s\n",
            (double)update_time);
    fprintf(stdout, "Load time: %.2f s\n",
            (double)load_time);
    fprintf(stdout, "Coding time: %.2lf s\n",
            (double)decode_time);
    fprintf(stdout, "Write time: %.2f s\n",
            (double)write_time);
#endif
    fprintf(stdout, "Total decompression time: %.2f s\n",
            (double)dec_time);

    return res;
}

/* --------------------------------------------------------------------------
 * Main program entry.
 */
static void usage(int err) {
    FILE *fp = err ? stderr : stdout;

    fprintf(fp, "Enano v%d.%d Author Guillermo Dufort y Alvarez, 2019-2020\n",
            MAJOR_VERS, MINOR_VERS);
    fprintf(fp, "The methods used for encoding the reads identifiers, and to model frequency counters, \n");
    fprintf(fp, "are the ones proposed by James Bonefield in FQZComp, with some modifications.\n");
    fprintf(fp, "The range coder is derived from Eugene Shelwien.\n\n");

    fprintf(fp, "To compress:\n  enano [options] [input_file [output_file]]\n\n");
    fprintf(fp, "    -c             To use MAX COMPRESION MODE. Default is FAST MODE.\n\n");
    fprintf(fp, "    -k <length>    Base sequence context length. Default is 7 (max 13).\n\n");
    fprintf(fp, "    -l <lenght>    Length of the DNA sequence context. Default is 6.\n\n");
    fprintf(fp, "    -t <num>       Maximum number of threads allowed to use by the compressor. Default is 8.\n\n");

    fprintf(fp, "To decompress:\n   enano -d [options] foo.enano foo.fastq\n");
    fprintf(fp, "    -t <num>       Maximum number of threads allowed to use by the decompressor. Default is 8.\n\n");

    exit(err);
}

int main (int argc, char **argv) {
    int res = 0;
    int decompress = 0;
    int opt;
    int in_fd = 0;
    int out_fd = 1;

    enano_params p;
    /* Initialise and parse command line arguments */
    p.klevel = DEFAULT_K_LEVEL;
    p.llevel = DEFAULT_L_LEVEL;
    p.num_threads = DEFAULT_THREADS_NUM;
    p.blk_upd_freq = DEFAULT_BLK_UPD_FREQ;
    p.blk_upd_thresh = DEFAULT_BLK_UPD_THRESH;
    p.max_comp = false;

    while ((opt = getopt(argc, argv, "hdk:l:t:cb:")) != -1) {
        switch (opt) {
            case 'h':
                usage(0);

            case 'd':
                decompress = 1;
                break;

            case 'k': {
                char *end;
                p.klevel = strtol(optarg, &end, 10);
                if (p.klevel < 0 || p.klevel > 13)
                    usage(1);
                break;
            }

            case 'l':
                p.llevel = atoi(optarg);
                break;

            case 't':
                p.num_threads = atoi(optarg);
                break;

            case 'b':
                p.blk_upd_thresh = atoi(optarg);
                break;

            case 'c':
                p.max_comp = true;
                break;

            default:
                usage(1);
        }
    }

    if (argc - optind > 2 || argc <= 2 || (decompress == 1 && argc <= 3))
        usage(1);

    if (optind != argc) {
        int open_flag = O_RDONLY;
        if ((in_fd = open(argv[optind], open_flag)) == -1) {
            perror(argv[optind]);
            exit(1);
        }
        optind++;
    }

    if (optind != argc) {
        out_fd = open(argv[optind], O_RDWR | O_CREAT | O_TRUNC, 0666);
        if (out_fd == -1) {
            perror(argv[optind]);
            exit(1);
        }
        optind++;
    }

    in_buf = new char[BLK_SIZE];
    omp_set_num_threads(p.num_threads);

    if (decompress) {
#ifdef __DEBUG_LOG__
        fp_log_debug = fopen("decode_log.txt", "wt");
#endif
        unsigned char magic[9];

        /* Check magic number */
        if (9 != read(in_fd, magic, 9)) {
            fprintf(stderr, "Abort: truncated read.\n");
            return 1;
        }
        if (memcmp(".ena", magic, 4) != 0) {
            fprintf(stderr, "Unrecognised file format.\n");
            return 1;
        }
        if (magic[4] != MAJOR_VERS) {
            fprintf(stderr, "Unsupported file format version %d.%d\n",
                    magic[4], magic[5]);
            return 1;
        }

        p.klevel = magic[5] & 0x0f;
        if (p.klevel > 13 || p.klevel < 1) {
            fprintf(stderr, "Unexpected quality compression level %d\n",
                    p.klevel);
            return 1;
        }

        p.llevel = magic[6] & 0x0f;

        p.max_comp = magic[7] & 0x0f;

        p.blk_upd_thresh = magic[8] & 0xff;

        printf("Parameters - k: %d, l: %d, b: %d \n", p.klevel, p.llevel, p.blk_upd_thresh);

        B_CTX_LEN = p.llevel;
        B_CTX = (1 << (B_CTX_LEN * A_LOG));
        AVG_CANT = (B_CTX * Q_CTX);
        NS_MODEL_SIZE = pow5[p.klevel];

        if (p.max_comp)
            res = decode_st(in_fd, out_fd, &p);
        else
            res = decode(in_fd, out_fd, &p);
#ifdef __DEBUG_LOG__
        fclose(fp_log_debug);
#endif
    } else {
#ifdef __DEBUG_LOG__
        fp_log_debug = fopen("encode_log.txt", "wt");
#endif
        unsigned char magic[9] = {'.', 'e', 'n', 'a',
                                  MAJOR_VERS,
                                  (unsigned char) p.klevel, (unsigned char) p.llevel, (unsigned char) p.max_comp, (unsigned char) p.blk_upd_thresh
        };

        if (9 != write(out_fd, magic, 9)) {
            fprintf(stderr, "Abort: truncated write.\n");
            return 1;
        }

        printf("Parameters - k: %d, l: %d, b: %d \n", p.klevel, p.llevel, p.blk_upd_thresh);

        B_CTX_LEN = p.llevel;
        B_CTX = (1 << (B_CTX_LEN * A_LOG));
        AVG_CANT = (B_CTX * Q_CTX);
        NS_MODEL_SIZE = pow5[p.klevel];

        if (p.max_comp)
            res = encode_st(in_fd, out_fd, &p);
        else
            res = encode(in_fd, out_fd, &p);

#ifdef __DEBUG_LOG__
        fclose(fp_log_debug);
#endif
    }

    delete in_buf;
    return res;
}
