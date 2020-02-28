//
// Created by gdufort on 6/28/19.
//

#define BLK_SIZE 10000000 //10 MB

//Default parameters
#define DEFAULT_K_LEVEL 7
#define DEFAULT_L_LEVEL 6
#define DEFAULT_THREADS_NUM 8
#define DEFAULT_BLK_UPD_THRESH 32
#define DEFAULT_BLK_UPD_FREQ 4

//#define __TIMING__
//#define __ORDER_SYMBOLS__

#define MAJOR_VERS 1
#define MINOR_VERS 0

#define SEQ_ALPHA_SIZE 5
static const unsigned int pow5[] = {
        1,
        5,
        25,
        125,
        625,
        3125,
        15625,
        78125,
        390625,
        1953125,
        9765625,
        48828125,
        244140625,
        1220703125
};

#define DUPLICATE_NAME_LINES
//#define __GLOBAL_STATS__
//#define __DEBUG_LOG__
//#define __CONTEXT_STATS__

#define ROUND_UP(x)    (x == 0? 0 : 1 << ((x) - 1))
/* Quantize quality scores q' = ((q + ROUND_UP(QUANT_Q)) >> QUANT_Q) */

#define DIV_ROUND(q, SHIFT) ((q + ROUND_UP(SHIFT)) >> SHIFT)

/* Keep as a power of 2 */
#define QMAX 128
#define A_LOG 2

#define QUANT_Q 3

#define LOG_AVGS 4
#define LOG_AVGS_ERRS 2

#define DIF_CTX_LOG 3
#define DIF_CANT (1 << DIF_CTX_LOG)

#define QUANT_D_LOG 5
#define QUANT_D_CANT (1 << QUANT_D_LOG)
#define QUANT_D_MASK (QUANT_D_CANT - 1)
#define QUANT_D_MAX QUANT_D_MASK

#define Q_LOG (7 - QUANT_Q)
#define TOTAL_Q_LOG (DIF_CTX_LOG + Q_LOG)

#define Q_CTX (1 << TOTAL_Q_LOG)

#define BC_CTX (1 << (LOG_AVGS + LOG_AVGS_ERRS))

#define CTX_CNT (BC_CTX * Q_CTX)

#define Q_LOG_CANT (1 << Q_LOG)

#define AVG_SHIFT 4
#define TOTAL_ERR_SHIFT 4

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define ABS(N) ((N<0)?(-N):(N))
