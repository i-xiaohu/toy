//
// Created by ixiaohu on 2019/9/19.
//

#ifndef TOY_SAMOP_H
#define TOY_SAMOP_H

#include "kvec.h"

// SAM flag
#define SAMF_MUL_SEG    0x1   // template having multiple segments in sequencing
static inline int mul_seg(int flag) { return ((flag & SAMF_MUL_SEG) != 0); }
#define SAMF_BOTH_ALI   0x2   // each segment properly aligned according to the aligner
static inline int both_ali(int flag) { return ((flag & SAMF_BOTH_ALI) != 0); }
#define SAMF_UNMAP      0x4   // segment unmapped
static inline int unmap(int flag) { return ((flag & SAMF_UNMAP) != 0); }
#define SAMF_NEXT_UNMAP 0x8   // next segment in the template unmapped
static inline int next_unmap(int flag) { return ((flag & SAMF_NEXT_UNMAP) != 0); }
#define SAMF_RC         0x10  // SEQ being reverse complemented
static inline int is_rc(int flag) { return ((flag & SAMF_RC) != 0); }
#define SAMF_NEXT_RC    0x20  // SEQ of the next segment in the template being reverse complemented
static inline int next_rc(int flag) { return ((flag & SAMF_NEXT_RC) != 0); }
#define SAMF_READ1      0x40  // the first segment in the template
static inline int is_read1(int flag) { return ((flag & SAMF_READ1) != 0); }
#define SAMF_READ2      0x80  // the last segment in the template
static inline int is_read2(int flag) { return ((flag & SAMF_READ2) != 0); }
#define SAMF_SEC_ALI    0x100 // secondary alignment, multiple mapping
static inline int sec_ali(int flag) { return ((flag & SAMF_SEC_ALI) != 0); }
#define SAMF_N0_FLT     0x200 // not passing filters, such as platform/vendor quality controls
static inline int no_flt(int flag) { return ((flag & SAMF_N0_FLT) != 0); }
#define SAMF_PCR        0x400 // PCR or optical duplicate
static inline int pcr(int flag) { return ((flag & SAMF_PCR) != 0); }
#define SAMF_SUP_ALI    0x800 // supplementary alignment, chimeric alignment
static inline int sup_ali(int flag) { return ((flag & SAMF_SUP_ALI) != 0); }

typedef struct {
	int nm; // edit distance
	int as; // alignment score
} sam_opt_t;

// SAM record
typedef struct {
	char *qname;
	int flag;
	char *rname;
	int pos;
	int mapq;
	char *cigar;
	char *rnext;
	int pnext;  // position of next read
	int tlen;   // inferred insert size
	char *seq;
	char *qual;
	sam_opt_t opt;
	char *data; // All fields are allocated in one memory chunk.
	// todo: get options
} sam_core1_t;
typedef kvec_t(sam_core1_t) sam_core1_v;

// SAM header, only saves the information I focus on
typedef struct {
	char *VN;
	char *SO;
} HD_t;
typedef struct {
	char *SN;
	int LN;
} SQ_t;
typedef kvec_t(SQ_t) SQ_v;
typedef struct {
	char *ID;
	char *PN;
	char *VN;
	char *CL;
} PG_t;
typedef struct {
	HD_t hd;
	SQ_v sqv;
	PG_t pg;
} sam_hdr_t;

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * Extract the SAM-header from a SAM head line.
	 * @param line SAM head line.
	 * @param len The length of $line.
	 * @param h Returned header.
	 */
	void sam_header(const char *line, int len, sam_hdr_t *h);

	/**
	 * Extract the SAM-record from a SAM line.
	 * @param line Input SAM line.
	 * @param len The length of $line.
	 * @param r Obtained record. Warn: $r uses outside memory of $line.
	 */
	void sam_record1(char *line, int len, sam_core1_t *r);

#ifdef __cplusplus
}
#endif

#endif //TOY_SAMOP_H
