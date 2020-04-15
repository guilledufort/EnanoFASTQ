//
// Created by Guillermo Dufort y Álvarez on 9/17/19.
//
#include "parameters.h"

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>

#include <inttypes.h>

#ifndef WIN32

#include <unistd.h>

#else
#include "unistd.h"
#endif

#include <ctype.h>
#include <assert.h>

#include <math.h>
#include <errno.h>
#include <time.h>

/* Range Coder:
 * This is using Eugene Shelwien's code from coders6c2.zip.
 */
#include "clr.h"

/*
 * Order 0 models, optimsed for various sizes of alphabet.
 * order0_coder is the original Dmitry Shkarin code, tweaked a bit to
 * make RangeCoder a parameter and to work as a template for adjusting
 * the alphabet size.
 *
 * Simple_model is James Bonefield implementation for FQZComp, adhering to the same API as
 * Dmitry Shkarin's original ORDER_0_CODER model.
 *
 * Base model N is James Bonefield's model specialising in symbols 0, 1, 2, 3, with
 * dedicated code for handling N (it encodes whichever of 0, 1, 2, 3
 * is most common and returns that value).
 */
#ifdef ORIG_MODEL
#    define ORDER_0_CODER SIMPLE_MODEL
#    include "order0_coder.h"
#else

#include "simple_model.h" // SIMPLE_MODEL

#endif

#ifdef FORCE_N_TO_Q0
#include "base_model.h"       // BASE_MODEL
#define SEQ_ALPHA_SIZE 4
#else

#include "base_modelN.h"       // BASE_MODEL

#endif


#define DECODE_INT(a) ((a)[0] + ((a)[1]<<8) + ((a)[2]<<16) + ((a)[3]<<24))
#define UPDATE_CONTEXT(ctx, b) ((ctx*5 + b) % NS_MODEL_SIZE)
/*
 * enano parameter block.
 */
typedef struct {
    uint8_t klevel;        // -s level
    uint8_t llevel;     // -l length
    uint8_t num_threads;
    bool max_comp;
    uint8_t blk_upd_freq;
    uint8_t blk_upd_thresh;
} enano_params;

typedef struct {
    //Read lengths
    SIMPLE_MODEL<256> model_len1;
    SIMPLE_MODEL<256> model_len2;
    SIMPLE_MODEL<256> model_len3;
    SIMPLE_MODEL<2> model_same_len;
    BASE_MODEL<uint8_t>* model_seq8;

    // Names
    SIMPLE_MODEL<256> model_name_prefix[256];
    SIMPLE_MODEL<256> model_name_suffix[256];
    SIMPLE_MODEL<256> model_name_len[256];
    SIMPLE_MODEL<128> model_name_middle[8192];

    // Qualities
    SIMPLE_MODEL<QUANT_D_CANT> model_qual_quant[CTX_CNT];
    SIMPLE_MODEL<QMAX - QUANT_D_CANT> quant_top;
} context_models;

/*
 * The enano class itself
 */
class Compressor {
public:
    Compressor();

    Compressor(enano_params *p);

    ~Compressor();

    void soft_reset();

    //Copies the statistical stats from cmp to this
    void copy_stats(context_models* ctx_m);

    /* Compression metrics */
    uint64_t base_in, base_out;
    uint64_t qual_in, qual_out;
    uint64_t name_in, name_out;
    uint64_t total_in, total_out;

    void compress_r1();

    void compress_r2();

    void compress_r3();

    void decompress_r1();

    void decompress_r2();

    void decompress_r3();

    bool output_block(int out_fd);
//protected:
    /* --- Parameters passed into the constructor */
    bool updateModel, maxCompression;

    context_models * cm;

    uint AVG_CANT, B_CTX, B_MASK, B_CTX_LEN, NS_MODEL_SIZE;

    char not_nl[256];

    char* decode_buf;

    int L[256];          // Sequence table lookups ACGTN->0..4

    unsigned char QDif[8][14] = {
            {0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7},
            {1, 0, 2, 3, 4, 5, 5, 6, 6, 6, 7},
            {2, 1, 0, 3, 4, 5, 5, 6, 6, 6, 7},
            {3, 2, 1, 0, 3, 4, 5, 5, 6, 6, 7},
            {3, 3, 2, 1, 0, 4, 5, 5, 6, 6, 7},
            {3, 3, 3, 2, 1, 0, 4, 5, 6, 6, 7},
            {3, 3, 3, 2, 2, 1, 0, 4, 5, 6, 7},
            {0, 1, 2, 3, 4, 3, 4, 3, 4, 5, 6}
    };

    unsigned char QBin[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14,
                                    14, 14, 14, 14, 14, 14, 14, 15,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,
                                    15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,
                                    15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,
                                    15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,
                                    15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15};

    /* --- Buffers */
    // Input and output buffers; Output buffer needs to be more than BLK_SIZE
    char out_buf[BLK_SIZE + BLK_SIZE / 2];
    int out_ind; // index into out_buf.
    int comp_len;
    int uncomp_len;

    int ns;
    int seq_len;
    char name_buf[BLK_SIZE/4];
    char seq_buf[BLK_SIZE];
    char qual_buf[BLK_SIZE];
    int name_len_a[BLK_SIZE / 9];
    int seq_len_a[BLK_SIZE / 9];
    char out0[BLK_SIZE]; // seq_len
    char out1[BLK_SIZE]; // name
    char out2[BLK_SIZE]; // seq
    char out3[BLK_SIZE]; // qual
    int sz0, sz1, sz2, sz3;
    char *in_buf0, *in_buf1, *in_buf2, *in_buf3;

    /* --- Models */
    // Sequence length
    int last_len;

    void encode_len(RangeCoder *rc, int len);

    int decode_len(RangeCoder *rc);

    char last_name[1024]; // Last name
    int last_name_len;    // Length of last name
    int last_p_len;       // Length of last common prefix
    int last_s_len;       // Length of last common suffix

    void encode_name(RangeCoder *rc, char *name, int len);

    int decode_name(RangeCoder *rc, char *name);

    void encode_seq8(RangeCoder *rc, char *seq, int len);

    void decode_seq8(RangeCoder *rc, char *seq, int len);

    // Quality
    uint16_t* ctx_avgs_sums;
    uint16_t* ctx_avgs_err_sums;
    uint32_t* ctx_err_avgs_total;

#ifdef __CONTEXT_STATS__
    uint ** context_stats;
    uint *quant_top_stats;
    uint *quant_bottom_stats;
#endif

    inline uint get_context(unsigned char s, unsigned char q1, unsigned char q2, uint &s_prev_ctx, uint &Q_prev_ctx);

    void encode_qual(RangeCoder *rc, char *seq, char *qual, int len);

    void decode_qual(RangeCoder *rc, char *seq, char *qual, int len);
    /* --- Main functions for compressing and decompressing blocks */
    /* Parses full reads into buffers from in. Returns remainder read length.*/
    int fq_parse_reads(char* in_buf, int in_len,
                       char **in_end, int *nseqs, int *remainder_length, int &end);

    int fq_compress();

    void fq_decompress();

    void update_AccFreqs(context_models* ctx_m, bool decode);

};
