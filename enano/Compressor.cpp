
#include "Compressor.h"

#ifdef __DEBUG_LOG__
FILE *fp_log_debug = NULL;
unsigned int debug_count = 0;
#endif

#ifdef __GLOBAL_STATS__
extern double under_T;
extern double over_T;
#endif

/* -------------------------------------------------------------------------
 * Constructors and destructors
 */


Compressor::Compressor() {
    Compressor(NULL);
}

Compressor::Compressor(enano_params *p) {

    cm = new context_models;

    B_CTX_LEN = p->llevel;
    B_CTX = (1 << (B_CTX_LEN * A_LOG));
    AVG_CANT = (B_CTX * Q_CTX);
    B_MASK = (B_CTX - 1);
    NS_MODEL_SIZE = pow5[p->klevel];

    for (uint i = 0; i < 256; i++)
        not_nl[i] = 1;
    not_nl['\r'] = not_nl['\n'] = 0;

    updateModel = true;
    maxCompression = p->max_comp;

    /* ACGTN* */
    for (int i = 0; i < 256; i++)
        L[i] = 0;

    L['A'] = L['a'] = 0;
    L['C'] = L['c'] = 1;
    L['G'] = L['g'] = 2;
    L['T'] = L['t'] = 3;
    L['N'] = L['n'] = 4;

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

    memset(ctx_err_avgs_total, 1 << TOTAL_ERR_SHIFT, Q_CTX * sizeof(uint16_t));

#ifdef __CONTEXT_STATS__
    context_stats = new uint *[CTX_CNT];

    for (int i = 0; i < CTX_CNT; i++) {
        context_stats[i] = new uint[1 + QUANT_D_CANT];
        memset(context_stats[i], 0, (1 + QUANT_D_CANT) * sizeof(uint));
    }

#endif

    /* Name settings */
    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;

    name_in = name_out = 0;
    base_in = base_out = 0;
    qual_in = qual_out = 0;
}

void Compressor::update_AccFreqs(context_models* ctx_m, bool decode){
    uint i;
    for (i = 0; i < 256; i++) {
        ctx_m->model_name_prefix[i].updateModelAccFrecs(decode);
        ctx_m->model_name_suffix[i].updateModelAccFrecs(decode);
        ctx_m->model_name_len[i].updateModelAccFrecs(decode);
    }

    for (i = 0; i < 8192; i++)
        ctx_m->model_name_middle[i].updateModelAccFrecs(decode);

    for (i = 0; i < CTX_CNT; i++)
        ctx_m->model_qual_quant[i].updateModelAccFrecs(decode);

    ctx_m->quant_top.updateModelAccFrecs(decode);
    ctx_m->model_len1.updateModelAccFrecs(decode);
    ctx_m->model_len2.updateModelAccFrecs(decode);
    ctx_m->model_len3.updateModelAccFrecs(decode);
    ctx_m->model_same_len.updateModelAccFrecs(decode);
}

void Compressor::soft_reset() {
    /* Name settings */
    memset(last_name, ' ', 1024);
    last_name_len = 0;
    last_p_len = 0;
    last_s_len = 0;

    /* Length settings */
    last_len = 0;
}

void Compressor::copy_stats(context_models* ctx_m){
    //Save pointer cause its going to get modified by the next line
    BASE_MODEL<uint8_t>* model_seq8_ptr = cm->model_seq8;
    *cm = *(ctx_m);
    //Copy memory
    memcpy(model_seq8_ptr, ctx_m->model_seq8, sizeof(BASE_MODEL<uint8_t>) * NS_MODEL_SIZE);
    //Set pointer again
    cm->model_seq8 = model_seq8_ptr;
}

Compressor::~Compressor() {
    delete [] ctx_avgs_sums;
    delete [] ctx_avgs_err_sums;
    delete [] ctx_err_avgs_total;
}

/* -------------------------------------------------------------------------
 * Name model
 */
void Compressor::encode_name(RangeCoder *rc, char *name, int len) {
    int p_len, s_len; // prefix and suffix length
    int i, j, k, last_char;

    // Prefix
    for (i = 0; i < len && i < last_name_len; i++) {
        if (name[i] != last_name[i])
            break;
    }
    p_len = i;

    // Suffix
    for (i = len - 1, j = last_name_len - 1; i >= 0 && j >= 0; i--, j--) {
        if (name[i] != last_name[j])
            break;
    }
    s_len = len - 1 - i;
    if (len - s_len - p_len < 0)
        s_len = len - p_len;

    if (updateModel) {
        if (maxCompression) {
            cm->model_name_prefix[last_p_len].encodeSymbolOrder(rc, p_len);
            cm->model_name_suffix[last_s_len].encodeSymbolOrder(rc, s_len);
            cm->model_name_len[last_name_len].encodeSymbolOrder(rc, len);
        } else {
            cm->model_name_prefix[last_p_len].encodeSymbol(rc, p_len);
            cm->model_name_suffix[last_s_len].encodeSymbol(rc, s_len);
            cm->model_name_len[last_name_len].encodeSymbol(rc, len);
        }
    } else {
        cm->model_name_prefix[last_p_len].encodeSymbolNoUpdate(rc, p_len);
        cm->model_name_suffix[last_s_len].encodeSymbolNoUpdate(rc, s_len);
        cm->model_name_len[last_name_len].encodeSymbolNoUpdate(rc, len);
    }

    last_p_len = p_len;
    last_s_len = s_len;

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;

        if (updateModel)
            if (maxCompression)
                cm->model_name_middle[last_char].encodeSymbolOrder(rc, name[i] & 0x7f);
            else
                cm->model_name_middle[last_char].encodeSymbol(rc, name[i] & 0x7f);
        else
            cm->model_name_middle[last_char].encodeSymbolNoUpdate(rc, name[i] & 0x7f);

        if (name[i] == ' ' && last_name[j] != ' ') j++;
        if (name[i] != ' ' && last_name[j] == ' ') j--;
        if (name[i] == ':' && last_name[j] != ':') j++;
        if (name[i] != ':' && last_name[j] == ':') j--;

        if (name[i] == ':' || name[i] == ' ') k = (k + 3) >> 2 << 2;

        lc2 = name[i] == last_name[j];
    }

    memcpy(last_name, name, len);
    last_name_len = len;
}

int Compressor::decode_name(RangeCoder *rc, char *name) {
    int p_len, s_len, len; // prefix and suffix length
    int i, j, k;
    int last_char;

    if (updateModel) {
        if (maxCompression) {
            p_len = cm->model_name_prefix[last_p_len].decodeSymbolOrder(rc);
            s_len = cm->model_name_suffix[last_s_len].decodeSymbolOrder(rc);
            len = cm->model_name_len[last_name_len].decodeSymbolOrder(rc);
        } else {
            p_len = cm->model_name_prefix[last_p_len].decodeSymbol(rc);
            s_len = cm->model_name_suffix[last_s_len].decodeSymbol(rc);
            len = cm->model_name_len[last_name_len].decodeSymbol(rc);
        }
    } else {
        p_len = cm->model_name_prefix[last_p_len].decodeSymbolNoUpdate(rc);
        s_len = cm->model_name_suffix[last_s_len].decodeSymbolNoUpdate(rc);
        len = cm->model_name_len[last_name_len].decodeSymbolNoUpdate(rc);
    }

    last_p_len = p_len;
    last_s_len = s_len;

    for (i = 0; i < p_len; i++)
        name[i] = last_name[i];

    int len2 = len - s_len, lc2 = p_len ? 1 : 0;
    for (i = j = p_len, k = 0; i < len2; i++, j++, k++) {
        unsigned char c;

        last_char = ((last_name[j] - 32) * 2 + lc2 + k * 64) % 8192;

        if (updateModel)
            if (maxCompression)
                c = cm->model_name_middle[last_char].decodeSymbolOrder(rc);
            else
                c = cm->model_name_middle[last_char].decodeSymbol(rc);
        else
            c = cm->model_name_middle[last_char].decodeSymbolNoUpdate(rc);

        //c = 'x';
        name[i] = c;

        if (c == ' ' && last_name[j] != ' ') j++;
        if (c != ' ' && last_name[j] == ' ') j--;
        if (c == ':' && last_name[j] != ':') j++;
        if (c != ':' && last_name[j] == ':') j--;

        if (name[i] == ':' || name[i] == ' ') k = (k + 3) >> 2 << 2;

        lc2 = c == last_name[j];
    }

    for (j = last_name_len - s_len; i < len; i++, j++)
        name[i] = last_name[j];

    memcpy(last_name, name, len);
    last_name_len = len;

    return len;
}

/* -------------------------------------------------------------------------
 * Sequence length model
 */
void Compressor::encode_len(RangeCoder *rc, int len) {
    if (updateModel) {
        if (maxCompression) {
            if (len != last_len) {
                cm->model_same_len.encodeSymbolOrder(rc, 0);
                cm->model_len1.encodeSymbolOrder(rc, len & 0xff);
                cm->model_len2.encodeSymbolOrder(rc, (len >> 8) & 0xff);
                cm->model_len3.encodeSymbolOrder(rc, (len >> 16) & 0xff);
            } else {
                cm->model_same_len.encodeSymbolOrder(rc, 1);
            }
        } else {
            if (len != last_len) {
                cm->model_same_len.encodeSymbol(rc, 0);
                cm->model_len1.encodeSymbol(rc, len & 0xff);
                cm->model_len2.encodeSymbol(rc, (len >> 8) & 0xff);
                cm->model_len3.encodeSymbol(rc, (len >> 16) & 0xff);
            } else {
                cm->model_same_len.encodeSymbol(rc, 1);
            }
        }
    } else {
        if (len != last_len) {
            cm->model_same_len.encodeSymbolNoUpdate(rc, 0);
            cm->model_len1.encodeSymbolNoUpdate(rc, len & 0xff);
            cm->model_len2.encodeSymbolNoUpdate(rc, (len >> 8) & 0xff);
            cm->model_len3.encodeSymbolNoUpdate(rc, (len >> 16) & 0xff);
        } else {
            cm->model_same_len.encodeSymbolNoUpdate(rc, 1);
        }
    }
}

int Compressor::decode_len(RangeCoder *rc) {
    if (updateModel) {
        if (maxCompression) {
            if (cm->model_same_len.decodeSymbolOrder(rc)) {
                return last_len;
            } else {
                int l1 = cm->model_len1.decodeSymbolOrder(rc);
                int l2 = cm->model_len2.decodeSymbolOrder(rc);
                int l3 = cm->model_len3.decodeSymbolOrder(rc);
                last_len = l1 + (l2 << 8) + (l3 << 16);
                return last_len;
            }
        } else {
            if (cm->model_same_len.decodeSymbol(rc)) {
                return last_len;
            } else {
                int l1 = cm->model_len1.decodeSymbol(rc);
                int l2 = cm->model_len2.decodeSymbol(rc);
                int l3 = cm->model_len3.decodeSymbol(rc);
                last_len = l1 + (l2 << 8) + (l3 << 16);
                return last_len;
            }
        }
    } else {
        if (cm->model_same_len.decodeSymbolNoUpdate(rc)) {
            return last_len;
        } else {
            int l1 = cm->model_len1.decodeSymbolNoUpdate(rc);
            int l2 = cm->model_len2.decodeSymbolNoUpdate(rc);
            int l3 = cm->model_len3.decodeSymbolNoUpdate(rc);
            last_len = l1 + (l2 << 8) + (l3 << 16);
            return last_len;
        }
    }
}


/* -------------------------------------------------------------------------
 * Sequence model
 */
void Compressor::encode_seq8(RangeCoder *rc, char *seq, int len) {
    int last;
    // Corresponds to a sequence of NS consecutive 'N'
    last = NS_MODEL_SIZE - 1;

    for (int i = 0; i < len; i++) {

        unsigned char b = L[(unsigned char) seq[i]];
        if (updateModel)
            cm->model_seq8[last].encodeSymbol(rc, b);
        else
            cm->model_seq8[last].encodeSymbolNoUpdate(rc, b);

        last = UPDATE_CONTEXT(last, b);
    }

}

void Compressor::decode_seq8(RangeCoder *rc, char *seq, int len) {
    int last;
    const char *dec = "ACGTN";

    // Corresponds to a sequence of NS consecutive 'N'
    last = NS_MODEL_SIZE - 1;

    for (int i = 0; i < len; i++) {
        unsigned char b;
        if (updateModel)
            b = cm->model_seq8[last].decodeSymbol(rc);
        else
            b = cm->model_seq8[last].decodeSymbolNoUpdate(rc);
        *seq++ = dec[b];
        last = UPDATE_CONTEXT(last, b);
    }

}

/* -------------------------------------------------------------------------
 * Quality model
 */

inline uint Compressor::get_context(unsigned char b, unsigned char q1, unsigned char q2, uint &B_prev_ctx, uint &Q_prev_ctx) {

    int err;
    uint nctx;

    //Update averages
    nctx = (B_prev_ctx << TOTAL_Q_LOG) + Q_prev_ctx;
    err = q1 - (DIV_ROUND(ctx_avgs_sums[nctx], AVG_SHIFT));
    ctx_avgs_sums[nctx] += err;
    uint abs_err = ABS(err);
    ctx_avgs_err_sums[nctx] += abs_err - (DIV_ROUND(ctx_avgs_err_sums[nctx], AVG_SHIFT));
    ctx_err_avgs_total[Q_prev_ctx] += abs_err - (DIV_ROUND(ctx_err_avgs_total[Q_prev_ctx], TOTAL_ERR_SHIFT));

    //We don't consider N for context. That saves memory.
    if (b == 'N')
        b = 'A';
    B_prev_ctx = ((B_prev_ctx << A_LOG) + L[b]) & B_MASK;

    //Current q_ctx
    uint q1_quant = QBin[q1];

    int dif_ctx = 0;
    if (q1 < 7) {
        dif_ctx = QDif[q1][MIN(q2, 10)];
    } else {
        int dif = q2 - q1;

        if (dif == 0) {
            dif_ctx = 0;
        } else {
            if (dif < 0)
                dif_ctx = QDif[7][MIN(-2 * dif - 1, 9)];
            else
                dif_ctx = QDif[7][MIN(2 * dif, 10)];
        }
    }

    Q_prev_ctx = ((dif_ctx << Q_LOG) + q1_quant);

    nctx = (B_prev_ctx << TOTAL_Q_LOG) + Q_prev_ctx;

    uint avg = QBin[DIV_ROUND(ctx_avgs_sums[nctx], AVG_SHIFT)];
    uint total_err_avg = (ctx_err_avgs_total[Q_prev_ctx] >> (TOTAL_ERR_SHIFT - AVG_SHIFT));
    err = ctx_avgs_err_sums[nctx];
    uint err_c = 0;

    if (err < (total_err_avg >> 1))
        err_c = 0;
    else if (err < total_err_avg)
        err_c = 1;
    else if (err <  (total_err_avg << 1))
        err_c = 2;
    else
        err_c = 3;

    return (err_c << (TOTAL_Q_LOG + LOG_AVGS)) + (avg << (TOTAL_Q_LOG)) + Q_prev_ctx;
}


void Compressor::encode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    int i, next_b;
    next_b = 1 + B_CTX_LEN/2;
    uint B_prev_ctx = 0, Q_prev_ctx = 0, ctx = 0;
    uint q1 = 0, q2 = 0;

    // Get first context
    for (i = 0; i < next_b; i++) {
        if (i < len)
            ctx = get_context(seq[i], q1, q2, B_prev_ctx, Q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, B_prev_ctx, Q_prev_ctx);
    }

    for (i = 0; i < len; i++, next_b++) {
        q1 = (qual[i] - '!') & (QMAX - 1);

        if (updateModel) {
            if (maxCompression) {
                if (q1 < QUANT_D_MAX) {
                    #ifdef __GLOBAL_STATS__
                    under_T++;
                    #endif
                    cm->model_qual_quant[ctx].encodeSymbolOrder(rc, q1);
                } else {
                    #ifdef __GLOBAL_STATS__
                    over_T++;
                    #endif
                    cm->model_qual_quant[ctx].encodeSymbolOrder(rc, QUANT_D_MAX);
                    cm->quant_top.encodeSymbolOrder(rc, q1 - QUANT_D_MAX);
                }
            } else {
                if (q1 < QUANT_D_MAX) {
                    #ifdef __GLOBAL_STATS__
                    under_T++;
                    #endif
                    cm->model_qual_quant[ctx].encodeSymbol(rc, q1);
                } else {
                    #ifdef __GLOBAL_STATS__
                    over_T++;
                    #endif
                    cm->model_qual_quant[ctx].encodeSymbol(rc, QUANT_D_MAX);
                    cm->quant_top.encodeSymbol(rc, q1 - QUANT_D_MAX);
                }
            }
        } else {
            if (q1 < QUANT_D_MAX) {
                #ifdef __GLOBAL_STATS__
                under_T++;
                #endif
                cm->model_qual_quant[ctx].encodeSymbolNoUpdate(rc, q1);
            } else {
                    #ifdef __GLOBAL_STATS__
                    over_T++;
                    #endif
                cm->model_qual_quant[ctx].encodeSymbolNoUpdate(rc, QUANT_D_MAX);
                cm->quant_top.encodeSymbolNoUpdate(rc, q1 - QUANT_D_MAX);
            }
        }
#ifdef  __CONTEXT_STATS__
        context_stats[ctx][0] += 1; //Add one to the total
        context_stats[ctx][MIN(q1+1, QUANT_D_MAX+1)] += 1;
#endif
#ifdef __DEBUG_LOG__
        fprintf(fp_log_debug, "ctx: %d\tq2: %d, q1: %d, dif_qual: %d\n", (int) ctx, (int) q2, (int) q1, (int) dif_qual);
        debug_count++;
#endif
        if (next_b < len)
            ctx = get_context(seq[next_b], q1, q2, B_prev_ctx, Q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, B_prev_ctx, Q_prev_ctx);

        q2 = q1;
    }
}

void Compressor::decode_qual(RangeCoder *rc, char *seq, char *qual, int len) {
    int i;
    int q1 = 0, q2 = 0;
    int next_b = 1 + B_CTX_LEN / 2;
    uint b_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;

    for (i = 0; i < next_b; i++) {
        if (i < len)
            ctx = get_context(seq[i], q1, q2, b_prev_ctx, q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, b_prev_ctx, q_prev_ctx);
    }

    for (i = 0; i < len; i++, next_b++) {
        unsigned q1;

        if (updateModel) {
            if (maxCompression) {
                q1 = (unsigned char) cm->model_qual_quant[ctx].decodeSymbolOrder(rc);
                if (q1 == QUANT_D_MAX) {
                    q1 += (unsigned char) cm->quant_top.decodeSymbolOrder(rc);
                }
            } else {
                q1 = (unsigned char) cm->model_qual_quant[ctx].decodeSymbol(rc);
                if (q1 == QUANT_D_MAX) {
                    q1 += (unsigned char) cm->quant_top.decodeSymbol(rc);
                }
            }
        } else {
            q1 = (unsigned char) cm->model_qual_quant[ctx].decodeSymbolNoUpdate(rc);
            if (q1 == QUANT_D_MAX) {
                q1 += (unsigned char) cm->quant_top.decodeSymbolNoUpdate(rc);
            }
        }
#ifdef __DEBUG_LOG__
        fprintf(fp_log_debug, "ctx: %d\tq2: %d, q1: %d, dif_qual: %d\n", (int) ctx, (int) q2, (int) q1, (int) dif_qual);
//        printf("cnt: %d\n", debug_count);//10050553
        debug_count++;
#endif

        if (next_b < len)
            ctx = get_context(seq[next_b], q1, q2, b_prev_ctx, q_prev_ctx);
        else
            ctx = get_context('A', q1, q2, b_prev_ctx, q_prev_ctx);

        qual[i] = q1 + '!';

        q2 = q1;
    }
}


/* --------------------------------------------------------------------------
 * Compression functions.
 */

/* Sequence length & name */
void Compressor::compress_r1() {
    char *name_p = name_buf;
    RangeCoder rc;

    rc.output(out1);
    rc.StartEncode();

    for (int i = 0; i < ns; i++) {
        encode_name(&rc, name_p, name_len_a[i]);
        name_p += name_len_a[i];
    }

    rc.FinishEncode();

    sz1 = rc.size_out();
    name_in += name_p - name_buf;
    name_out += sz1;
}

/* Sequence itself */
void Compressor::compress_r2() {
    char *seq_p = seq_buf;
    RangeCoder rc;

    rc.output(out2);
    rc.StartEncode();
    for (int i = 0; i < ns; i++) {
        encode_seq8(&rc, seq_p, seq_len_a[i]);
        seq_p += seq_len_a[i];
    }
    rc.FinishEncode();

    sz2 = rc.size_out();
    base_in += seq_p - seq_buf;
    base_out += sz2;
}

/* Quality values */
void Compressor::compress_r3() {

    char *qual_p = qual_buf;
    char *seq_p = seq_buf;
    RangeCoder rc;

    rc.output(out3);
    rc.StartEncode();

    for (int i = 0; i < ns; i++) {
        encode_qual(&rc, seq_p, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
        seq_p += seq_len_a[i];
    }

    rc.FinishEncode();

    sz3 = rc.size_out();
    qual_in += qual_p - qual_buf;
    qual_out += sz3;
}


/*
 * Reads from in[0..in_len-1] and writes a compressed copy to out, setting
 * *out_len to the returned size. The caller needs to ensure out is large
 * enough.
 *
 * We can only compress entire sequences and in[] may end in a partial
 * sequence. Hence *in_end is set to a pointer to the end of the processed
 * buffer.
 *
 * Also returns *nseqs filled out.
 *
 * Returns total compressed length on success, with data in out[]
 *        -1 on failure
 */
int Compressor::fq_parse_reads(char* in_buf, int in_len,
                               char **in_end, int *nseqs, int *remainder_length, int &end) {
    char *in = &in_buf[0];
    end = 0;
    int end_hash = 0;
    int i, j;

    char *name_p = name_buf;
    char *seq_p = seq_buf;
    char *qual_p = qual_buf;

    ns = 0;

    /* Parse and separate into name, seq, qual buffers */
    seq_len = 0;

    uint len;

    for (i = 0; i < in_len;) {
        /* Name */
        if (in[i] != '@')
            return -1;

        j = i;
        i++;
        while (i < in_len && not_nl[(uc) in[i]])
            i++;

        len = i - j - 1;
        name_len_a[ns] = len;
        memcpy(name_p, &in[j+1],len * sizeof(char));
        name_p += len;

        if (++i >= in_len)
            break;

        /* Sequence */
        for (j = i; i < in_len && not_nl[(uc) in[i]]; i++);

        len = i - j;
        seq_len_a[ns] = len;
        memcpy(seq_p, &in[j],len * sizeof(char));
        seq_p += len;

        if (++i >= in_len)
            break;

        /* +name, assume to be identical to @name */
        if (in[i] != '+')
            return -2;

        for (; i < in_len && not_nl[(uc) in[i]]; i++);
        if (++i >= in_len)
            break;

        /* Quality */
        if (i + seq_len_a[ns] > in_len){
            i = in_len + 1;
            break;
        }

        memcpy(qual_p, &in[i], seq_len_a[ns] * sizeof(char));
        qual_p += seq_len_a[ns];
        i += seq_len_a[ns];

        if (++i > in_len)
            break;

        end = i;
        end_hash = i;

        if (seq_len == 0)
            seq_len = seq_len_a[ns];
        else if (seq_len != seq_len_a[ns])
            seq_len = -1;

        ns++;
    }

    *in_end = in + end_hash;

    /* rl = first_pos - last_pos + 1
     (we dont add 1 because last pos is already one more than needed) */
    *remainder_length = i - end_hash;

    if (i == BLK_SIZE)
        *remainder_length += 1;

    *nseqs = ns;

    return 0;
}

int Compressor::fq_compress(){
    /* Encode seq len, we have a dependency on this for seq/qual */
    char *out = out_buf + 4;
    RangeCoder rc;
    rc.output(out0);
    rc.StartEncode();
    for (int i = 0; i < ns; i++) {
        encode_len(&rc, seq_len_a[i]);
    }
    rc.FinishEncode();
    sz0 = rc.size_out();

#pragma omp parallel sections
    {
#pragma omp section
        {
            compress_r1();
        }
#pragma omp section
        {
            compress_r2();
        }
#pragma omp section
        {
            compress_r3();
        }
    }

    /* Concatenate compressed output into a single block */
    char *out_p = out;

    *out_p++ = (ns >> 0) & 0xff;  /* Number of sequences */
    *out_p++ = (ns >> 8) & 0xff;
    *out_p++ = (ns >> 16) & 0xff;
    *out_p++ = (ns >> 24) & 0xff;

    *out_p++ = (sz0 >> 0) & 0xff;  /* Size of 4 range-coder blocks */
    *out_p++ = (sz0 >> 8) & 0xff;
    *out_p++ = (sz0 >> 16) & 0xff;
    *out_p++ = (sz0 >> 24) & 0xff;

    *out_p++ = (sz1 >> 0) & 0xff;
    *out_p++ = (sz1 >> 8) & 0xff;
    *out_p++ = (sz1 >> 16) & 0xff;
    *out_p++ = (sz1 >> 24) & 0xff;

    *out_p++ = (sz2 >> 0) & 0xff;
    *out_p++ = (sz2 >> 8) & 0xff;
    *out_p++ = (sz2 >> 16) & 0xff;
    *out_p++ = (sz2 >> 24) & 0xff;

    *out_p++ = (sz3 >> 0) & 0xff;
    *out_p++ = (sz3 >> 8) & 0xff;
    *out_p++ = (sz3 >> 16) & 0xff;
    *out_p++ = (sz3 >> 24) & 0xff;

#ifdef __DEBUG_BLOCKS__
    printf("ns : %ld\n", ns);
    printf("lens : %ld\n", sz0);
    printf("names : %ld\n", sz1);
    printf("bases : %ld\n", sz2);
    printf("quals : %ld\n", sz3);
#endif

    memcpy(out_p, out0, sz0);
    out_p += sz0;
    memcpy(out_p, out1, sz1);
    out_p += sz1;
    memcpy(out_p, out2, sz2);
    out_p += sz2;
    memcpy(out_p, out3, sz3);
    out_p += sz3;

    comp_len = out_p - out;

    return 0;
}
/* --------------------------------------------------------------------------
 * Decompression functions.
 */

void Compressor::decompress_r1(void) {
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    char *name_p = name_buf;

    for (int i = 0; i < ns; i++) {
        *name_p++ = '@';
        name_p += decode_name(&rc, name_p);
        *name_p++ = '\n';
    }

    rc.FinishDecode();
}

void Compressor::decompress_r2(void) {
    RangeCoder rc;
    rc.input(in_buf2);
    rc.StartDecode();

    char *seq_p = seq_buf;
    for (int i = 0; i < ns; i++) {
        decode_seq8(&rc, seq_p, seq_len_a[i]);
        seq_p += seq_len_a[i];
    }
    rc.FinishDecode();
}

void Compressor::decompress_r3(void) {
    RangeCoder rc;
    rc.input(in_buf3);
    rc.StartDecode();

    char *seq_p = seq_buf;
    char *qual_p = qual_buf;
    for (int i = 0; i < ns; i++) {
        decode_qual(&rc, seq_p, qual_p, seq_len_a[i]);
        qual_p += seq_len_a[i];
        seq_p += seq_len_a[i];
    }
    rc.FinishDecode();
}

/* Decompress a single block */
void Compressor::fq_decompress() {

    char *name_p, *seq_p, *qual_p;
    char *in = decode_buf;
    uint32_t nseqs = DECODE_INT((unsigned char *) (in));
    uint32_t sz0 = DECODE_INT((unsigned char *) (in + 4));
    uint32_t sz1 = DECODE_INT((unsigned char *) (in + 8));
    uint32_t sz2 = DECODE_INT((unsigned char *) (in + 12));
    uint32_t sz3 = DECODE_INT((unsigned char *) (in + 16));

#ifdef __DEBUG_BLOCKS__
    printf("ns : %ld\n", nseqs);
    printf("lens : %ld\n", sz0);
    printf("names : %ld\n", sz1);
    printf("bases : %ld\n", sz2);
    printf("quals : %ld\n", sz3);
#endif
    in += 20;
    ns = nseqs;

    in_buf0 = in;
    in += sz0;
    in_buf1 = in;
    in += sz1;
    in_buf2 = in;
    in += sz2;
    in_buf3 = in;
    in += sz3;

    RangeCoder rc0;
    rc0.input(in_buf0);
    rc0.StartDecode();

    for (int i = 0; i < ns; i++)
        seq_len_a[i] = decode_len(&rc0);
    rc0.FinishDecode();

#pragma omp parallel sections
    {
#pragma omp section
        {
            decompress_r1();
        }
#pragma omp section
        {
            decompress_r2();
            decompress_r3();
        }
    }

    /* Stick together the arrays into out_buf */
    out_ind = 0;
    name_p = name_buf;
    seq_p = seq_buf;
    qual_p = qual_buf;

    for (int i = 0; i < ns; i++) {
        /* name */
        char *aux_name_p = name_p;
        while ((out_buf[out_ind++] = *name_p++) != '\n');

        /* seq */
        for (int j = 0; j < seq_len_a[i]; j++)
            out_buf[out_ind++] = *seq_p++;

        out_buf[out_ind++] = '\n';
        out_buf[out_ind++] = '+';
#ifdef DUPLICATE_NAME_LINES
        while ((out_buf[out_ind++] = *(++aux_name_p)) != '\n');
#else
        out_buf[out_ind++] = '\n';
#endif
        /* qual */
        for (int j = 0; j < seq_len_a[i]; j++) {
            out_buf[out_ind++] = *qual_p++;
        }
        out_buf[out_ind++] = '\n';
    }

    uncomp_len = out_ind;
}

bool Compressor::output_block(int out_fd){
    out_buf[0] = (comp_len >> 0) & 0xff;
    out_buf[1] = (comp_len >> 8) & 0xff;
    out_buf[2] = (comp_len >> 16) & 0xff;
    out_buf[3] = (comp_len >> 24) & 0xff;
    comp_len += 4;

    int sz = write(out_fd, out_buf, comp_len);
    return (sz == comp_len) ;
}
