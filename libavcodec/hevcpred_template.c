/*
 * HEVC video Decoder
 *
 * Copyright (C) 2012 Guillaume Martres
 *
 * This file is part of Libav.
 *
 * Libav is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Libav is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Libav; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "libavutil/pixdesc.h"
#include "bit_depth_template.c"
#include "hevcpred.h"
//#define USE_SSE
#ifdef USE_SSE
#include <emmintrin.h>
#include <x86intrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#endif

#ifdef USE_SSE
#define SSE_planar
#endif

#define POS(x, y) src[(x) + stride * (y)]

static void FUNCC(intra_pred)(HEVCContext *s, int x0, int y0, int log2_size, int c_idx, int entry)
{
#define MIN_TB_ADDR_ZS(x, y)                                            \
    s->pps->min_tb_addr_zs[(y) * s->sps->pic_width_in_min_tbs + (x)]

#define EXTEND_LEFT(ptr, length)                \
    for (i = 0; i < (length); i++)              \
        (ptr)[-(i+1)] = (ptr)[0];
#define EXTEND_RIGHT(ptr, length)               \
    for (i = 0; i < (length); i++)              \
        (ptr)[i+1] = (ptr)[0];
#define EXTEND_UP(ptr, length) EXTEND_LEFT(ptr, length)
#define EXTEND_DOWN(ptr, length) EXTEND_RIGHT(ptr, length)

    int i;
    int hshift = s->sps->hshift[c_idx];
    int vshift = s->sps->vshift[c_idx];
    int size = (1 << log2_size);
    int size_in_luma = size << hshift;
    int size_in_tbs = size_in_luma >> s->sps->log2_min_transform_block_size;
    int x = x0 >> hshift;
    int y = y0 >> vshift;
    int x_tb = x0 >> s->sps->log2_min_transform_block_size;
    int y_tb = y0 >> s->sps->log2_min_transform_block_size;
    int cur_tb_addr = MIN_TB_ADDR_ZS(x_tb, y_tb);

    ptrdiff_t stride = s->frame->linesize[c_idx] / sizeof(pixel);
    pixel *src = (pixel*)s->frame->data[c_idx] + x + y * stride;

    enum IntraPredMode mode = c_idx ? s->pu.intra_pred_mode_c[entry] :
                              s->tu[entry].cur_intra_pred_mode;

    pixel left_array[2*MAX_TB_SIZE+1], filtered_left_array[2*MAX_TB_SIZE+1];
    pixel top_array[2*MAX_TB_SIZE+1], filtered_top_array[2*MAX_TB_SIZE+1];
    pixel *left = left_array + 1;
    pixel *top = top_array + 1;
    pixel *filtered_left = filtered_left_array + 1;
    pixel *filtered_top = filtered_top_array + 1;

    int x0b = x0 & ((1 << s->sps->log2_ctb_size) - 1);
    int y0b = y0 & ((1 << s->sps->log2_ctb_size) - 1);

    int cand_up   = (s->ctb_up_flag[entry] || y0b);
    int cand_left   = (s->ctb_left_flag[entry] || x0b);
    
    int bottom_left_available = cand_left && (y_tb + size_in_tbs) < s->sps->pic_height_in_min_tbs &&
                                cur_tb_addr > MIN_TB_ADDR_ZS(x_tb - 1, y_tb + size_in_tbs);
    int left_available = cand_left;
    int top_left_available = cand_left && cand_up;
    int top_available = cand_up;
    //FIXME : top_right_available can be available even if cand_up is not 
    int top_right_available = cand_up && (x_tb + size_in_tbs) < (s->end_of_tiles_x[entry]>>s->sps->log2_min_transform_block_size) &&
                              cur_tb_addr > MIN_TB_ADDR_ZS(x_tb + size_in_tbs, y_tb - 1);

    int bottom_left_size = (FFMIN(y0 + 2*size_in_luma, s->sps->pic_height_in_luma_samples) -
                            (y0 + size_in_luma)) >> vshift;
    int top_right_size = (FFMIN(x0 + 2*size_in_luma, s->sps->pic_width_in_luma_samples) -
                          (x0 + size_in_luma)) >> hshift;
    if (s->pps->constrained_intra_pred_flag == 1) {
        int pic_width_in_min_pu  = s->sps->pic_width_in_min_cbs * 4;
        int size_pu     =  1 << s->sps->log2_min_pu_size;
        int x_pu        = x0 >> s->sps->log2_min_pu_size;
        int y_pu        = y0 >> s->sps->log2_min_pu_size;
        int x0_pu       = x0 & (size_pu - 1);
        int y0_pu       = y0 & (size_pu - 1);
        int x_left_pu   = x0_pu == 0 ? x_pu - 1 : x_pu;
        int x_right_pu  = x0_pu >= size_pu ? x_pu + 1 : x_pu;
        int y_top_pu    = y0_pu == 0 ? y_pu - 1 : y_pu;
        int y_bottom_pu = y0_pu >= size_pu ? y_pu + 1 : y_pu;
        if (bottom_left_available == 1)
            bottom_left_available = s->ref->tab_mvf[x_left_pu + y_bottom_pu * pic_width_in_min_pu].is_intra;
        if (left_available == 1)
            left_available = s->ref->tab_mvf[x_left_pu + y_pu * pic_width_in_min_pu].is_intra;
        if (top_left_available == 1)
            top_left_available = s->ref->tab_mvf[x_left_pu + y_top_pu * pic_width_in_min_pu].is_intra;
        if (top_available == 1)
            top_available = s->ref->tab_mvf[x_pu + y_top_pu * pic_width_in_min_pu].is_intra;
        if (top_right_available == 1)
            top_right_available = s->ref->tab_mvf[x_right_pu + y_top_pu * pic_width_in_min_pu].is_intra;
    }
    // Fill left and top with the available samples
    if (bottom_left_available) {
        for (i = 0; i < bottom_left_size; i++) {
            left[size + i] = POS(-1, size + i);
        }
        for (; i < size; i++) {
            left[size + i] = POS(-1, size + bottom_left_size - 1);
        }
    }
    if (left_available) {
        for (i = 0; i < size; i++)
            left[i] = POS(-1, i);
    }
    if (top_left_available)
        left[-1] = POS(-1, -1);
    if (top_available && top_right_available && top_right_size == size) {
        memcpy(&top[0], &POS(0, -1), size * sizeof(pixel));
        memcpy(&top[size], &POS(size, -1), top_right_size * sizeof(pixel));
    } else {
        if (top_available)
            memcpy(&top[0], &POS(0, -1), size * sizeof(pixel));
        if (top_right_available) {
            memcpy(&top[size], &POS(size, -1), top_right_size * sizeof(pixel));
            for (i = top_right_size; i < size; i++)
                top[size + i] = POS(size + top_right_size - 1, -1);
        }
    }
    // Infer the unavailable samples
    if (!bottom_left_available) {
        if (left_available) {
            EXTEND_DOWN(&left[size-1], size);
        } else if (top_left_available) {
            EXTEND_DOWN(&left[-1], 2*size);
            left_available = 1;
        } else if (top_available) {
            left[-1] = top[0];
            EXTEND_DOWN(&left[-1], 2*size);
            top_left_available = 1;
            left_available = 1;
        } else if (top_right_available) {
            EXTEND_LEFT(&top[size], size);
            left[-1] = top[0];
            EXTEND_DOWN(&left[-1], 2*size);
            top_available = 1;
            top_left_available = 1;
            left_available = 1;
        } else { // No samples available
            top[0] = left[-1] = (1 << (BIT_DEPTH - 1));
            EXTEND_RIGHT(&top[0], 2*size-1);
            EXTEND_DOWN(&left[-1], 2*size);
        }
    }

    if (!left_available) {
        EXTEND_UP(&left[size], size);
    }
    if (!top_left_available) {
        left[-1] = left[0];
    }
    if (!top_available) {
        top[0] = left[-1];
        EXTEND_RIGHT(&top[0], size-1);
    }
    if (!top_right_available)
        EXTEND_RIGHT(&top[size-1], size);

    top[-1] = left[-1];

#undef EXTEND_LEFT
#undef EXTEND_RIGHT
#undef EXTEND_UP
#undef EXTEND_DOWN
#undef MIN_TB_ADDR_ZS

    // Filtering process
    if (c_idx == 0 && mode != INTRA_DC && size != 4) {
        int intra_hor_ver_dist_thresh[] = { 7, 1, 0 };
        int min_dist_vert_hor = FFMIN(FFABS((int)mode-26), FFABS((int)mode-10));
        if (min_dist_vert_hor > intra_hor_ver_dist_thresh[log2_size-3]) {
            int thresold = 1 << (BIT_DEPTH - 5);
            if (s->sps->sps_strong_intra_smoothing_enable_flag && log2_size == 5 &&
                FFABS(top[-1] + top[63] - 2 * top[31]) < thresold &&
                FFABS(left[-1] + left[63] - 2 * left[31]) < thresold) {
                // We can't just overwrite values in top because it could be a pointer into src
                filtered_top[-1] = top[-1];
                filtered_top[63] = top[63];
                for (i = 0; i < 63; i++) {
                    filtered_top[i] = ((64 - (i + 1))*top[-1] + (i + 1) * top[63] + 32) >> 6;
                }
                for (i = 0; i < 63; i++) {
                    left[i] = ((64 - (i + 1))*left[-1] + (i + 1) * left[63] + 32) >> 6;
                }
                top = filtered_top;
            } else {
                filtered_left[2*size-1] = left[2*size-1];
                filtered_top[2*size-1]  = top[2*size-1];
                for (i = 2*size-2; i >= 0; i--) {
                    filtered_left[i] = (left[i+1] + 2*left[i] + left[i-1] + 2) >> 2;
                }
                filtered_top[-1] = filtered_left[-1] = (left[0] + 2*left[-1] + top[0] + 2) >> 2;
                for (i = 2*size-2; i >= 0; i--) {
                    filtered_top[i] = (top[i+1] + 2*top[i] + top[i-1] + 2) >> 2;
                }
                left = filtered_left;
                top = filtered_top;
            }
        }
    }

    switch(mode) {
    case INTRA_PLANAR:
        s->hpc.pred_planar[ (log2_size - 2) ]((uint8_t*)src, (uint8_t*)top, (uint8_t*)left, stride, log2_size);
        break;
    case INTRA_DC:
        s->hpc.pred_dc((uint8_t*)src, (uint8_t*)top, (uint8_t*)left, stride, log2_size, c_idx);
        break;
    default:
        s->hpc.pred_angular((uint8_t*)src, (uint8_t*)top, (uint8_t*)left, stride, log2_size, c_idx,
                            mode);
        break;
    }

}

#ifdef SSE_planar
static void FUNCC(pred_planar_0)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    __m128i ly, t0, tx, l0, add, c0,c1,c2,c3, ly1, tmp1, tmp2, mask;
    t0= _mm_set1_epi16(top[4]);
    l0= _mm_set1_epi16(left[4]);
    add= _mm_set1_epi16(4);

	ly= _mm_loadu_si128((__m128i*)left);			//get 16 values
	ly= _mm_unpacklo_epi8(ly,_mm_setzero_si128());	//drop to 8 values 16 bit

	tx= _mm_loadu_si128((__m128i*)top);				//get 16 values
	tx= _mm_unpacklo_epi8(tx,_mm_setzero_si128());	//drop to 8 values 16 bit
	tmp1= _mm_set_epi16(0,0,0,0,0,1,2,3);
	tmp2= _mm_set_epi16(0,0,0,0,4,3,2,1);
	mask= _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1);
    for (y = 0; y < 4; y++){

    	ly1= _mm_set1_epi16(_mm_extract_epi16(ly,y));

    	c0= _mm_mullo_epi16(tmp1,ly1);
    	c1= _mm_mullo_epi16(tmp2,t0);
    	c2= _mm_mullo_epi16(_mm_set1_epi16(3 - y),tx);
    	c3= _mm_mullo_epi16(_mm_set1_epi16(1+y),l0);

    	c0= _mm_add_epi16(c0,c1);
    	c2= _mm_add_epi16(c2,c3);
    	c2= _mm_add_epi16(c2,add);
    	c0= _mm_add_epi16(c0,c2);

    	c0= _mm_srli_epi16(c0,3);

    	c0= _mm_packus_epi16(c0,_mm_setzero_si128());

    	_mm_maskmoveu_si128(c0,mask,(__m128i*)(src+y*stride)); //store only 4 values

    }
}
static void FUNCC(pred_planar_1)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;

    __m128i ly, t0, tx, l0, add, c0,c1,c2,c3, ly1, tmp1, tmp2;
    t0= _mm_set1_epi16(top[8]);
    l0= _mm_set1_epi16(left[8]);
    add= _mm_set1_epi16(8);

    ly= _mm_loadu_si128((__m128i*)left);			//get 16 values
    ly= _mm_unpacklo_epi8(ly,_mm_setzero_si128());	//drop to 8 values 16 bit

    tx= _mm_loadu_si128((__m128i*)top);				//get 16 values
    tx= _mm_unpacklo_epi8(tx,_mm_setzero_si128());	//drop to 8 values 16 bit
    tmp1= _mm_set_epi16(0,1,2,3,4,5,6,7);
    tmp2= _mm_set_epi16(8,7,6,5,4,3,2,1);

    for (y = 0; y < 8; y++){

    	ly1= _mm_set1_epi16(_mm_extract_epi16(ly,0));

    	c0= _mm_mullo_epi16(tmp1,ly1);
    	c1= _mm_mullo_epi16(tmp2,t0);
    	c2= _mm_mullo_epi16(_mm_set1_epi16(7 - y),tx);
    	c3= _mm_mullo_epi16(_mm_set1_epi16(1+y),l0);

    	c0= _mm_add_epi16(c0,c1);
    	c2= _mm_add_epi16(c2,c3);
    	c2= _mm_add_epi16(c2,add);
    	c0= _mm_add_epi16(c0,c2);

    	c0= _mm_srli_epi16(c0,4);

    	c0= _mm_packus_epi16(c0,_mm_setzero_si128());

    	_mm_storel_epi64((__m128i*)(src + y*stride), c0);	//store only 8

    	ly= _mm_srli_si128(ly,2);

    }

}
static void FUNCC(pred_planar_2)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;

    __m128i ly, t0, tx, l0, add, c0,c1,c2,c3, ly1, tmp1, tmp2,C0,C1,C2,C3;
    t0= _mm_set1_epi16(top[16]);
    l0= _mm_set1_epi16(left[16]);
    add= _mm_set1_epi16(16);

    ly= _mm_loadu_si128((__m128i*)left);			//get 16 values

    tx= _mm_loadu_si128((__m128i*)top);				//get 16 values
    tmp1= _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    tmp2= _mm_set_epi8(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);

    for (y = 0; y < 16; y++){
    	ly1= _mm_set1_epi16(_mm_extract_epi8(ly,0));

    	c0= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp1,_mm_setzero_si128()),ly1);
    	C0= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp1,_mm_setzero_si128()),ly1);

    	c1= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp2,_mm_setzero_si128()),t0);
    	C1= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp2,_mm_setzero_si128()),t0);

    	c2= _mm_mullo_epi16(_mm_set1_epi16(15 - y),_mm_unpacklo_epi8(tx,_mm_setzero_si128()));
    	C2= _mm_mullo_epi16(_mm_set1_epi16(15 - y),_mm_unpackhi_epi8(tx,_mm_setzero_si128()));

    	c3= _mm_mullo_epi16(_mm_set1_epi16(1+y),l0);

    	c0= _mm_add_epi16(c0,c1);
    	c2= _mm_add_epi16(c2,c3);
    	c2= _mm_add_epi16(c2,add);
    	c0= _mm_add_epi16(c0,c2);

    	C0= _mm_add_epi16(C0,C1);
    	C2= _mm_add_epi16(C2,c3);
    	C2= _mm_add_epi16(C2,add);
    	C0= _mm_add_epi16(C0,C2);

    	c0= _mm_srli_epi16(c0,5);
    	C0= _mm_srli_epi16(C0,5);

    	c0= _mm_packus_epi16(c0,C0);

    	_mm_storeu_si128((__m128i*)(src + y*stride), c0);
    	ly= _mm_srli_si128(ly,1);
    }

}
static void FUNCC(pred_planar_3)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;

    __m128i ly, LY, t0, tx, TX, l0, add, c0,c1,c2,c3, ly1, tmp1, tmp2, TMP1, TMP2,C0,C1,C2,C3;
        t0= _mm_set1_epi16(top[32]);
        l0= _mm_set1_epi16(left[32]);
        add= _mm_set1_epi16(32);

        ly= _mm_loadu_si128((__m128i*)left);			//get 16 values
        LY= _mm_loadu_si128((__m128i*)(left+16));

        tx= _mm_loadu_si128((__m128i*)top);				//get 16 values
        TX= _mm_loadu_si128((__m128i*)(top +16));
        TMP1= _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
        tmp1= _mm_set_epi8(16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
        tmp2= _mm_set_epi8(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);
        TMP2= _mm_set_epi8(32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17);

        for (y = 0; y < 16; y++){
        	//first half of 32
        	ly1= _mm_set1_epi16(_mm_extract_epi8(ly,0));
        	c3= _mm_mullo_epi16(_mm_set1_epi16(1+y),l0);

        	c0= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp1,_mm_setzero_si128()),ly1);
        	//printf("values check : tmp1 = %d, ly1= ");
        	C0= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp1,_mm_setzero_si128()),ly1);

        	c1= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp2,_mm_setzero_si128()),t0);
        	C1= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp2,_mm_setzero_si128()),t0);

        	c2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpacklo_epi8(tx,_mm_setzero_si128()));
        	C2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpackhi_epi8(tx,_mm_setzero_si128()));

        	c0= _mm_add_epi16(c0,c1);
        	c2= _mm_add_epi16(c2,c3);
        	c2= _mm_add_epi16(c2,add);
        	c0= _mm_add_epi16(c0,c2);

        	C0= _mm_add_epi16(C0,C1);
        	C2= _mm_add_epi16(C2,c3);
        	C2= _mm_add_epi16(C2,add);
        	C0= _mm_add_epi16(C0,C2);

        	c0= _mm_srli_epi16(c0,6);
        	C0= _mm_srli_epi16(C0,6);

        	c0= _mm_packus_epi16(c0,C0);

        	_mm_storeu_si128((__m128i*)(src + y*stride), c0);

        	// second half of 32

        	c0= _mm_mullo_epi16(_mm_unpacklo_epi8(TMP1,_mm_setzero_si128()),ly1);
        	C0= _mm_mullo_epi16(_mm_unpackhi_epi8(TMP1,_mm_setzero_si128()),ly1);

        	c1= _mm_mullo_epi16(_mm_unpacklo_epi8(TMP2,_mm_setzero_si128()),t0);
        	C1= _mm_mullo_epi16(_mm_unpackhi_epi8(TMP2,_mm_setzero_si128()),t0);

        	c2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpacklo_epi8(TX,_mm_setzero_si128()));
        	C2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpackhi_epi8(TX,_mm_setzero_si128()));



        	c0= _mm_add_epi16(c0,c1);
        	c2= _mm_add_epi16(c2,c3);
        	c2= _mm_add_epi16(c2,add);
        	c0= _mm_add_epi16(c0,c2);

        	C0= _mm_add_epi16(C0,C1);
        	C2= _mm_add_epi16(C2,c3);
        	C2= _mm_add_epi16(C2,add);
        	C0= _mm_add_epi16(C0,C2);

        	c0= _mm_srli_epi16(c0,6);
        	C0= _mm_srli_epi16(C0,6);

        	c0= _mm_packus_epi16(c0,C0);

        	_mm_storeu_si128((__m128i*)(src + 16 + y*stride), c0);

        	ly= _mm_srli_si128(ly,1);
        }

        for (y = 16; y < 32; y++){
        	//first half of 32
        	ly1= _mm_set1_epi16(_mm_extract_epi8(LY,0));
        	c3= _mm_mullo_epi16(_mm_set1_epi16(1+y),l0);

        	c0= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp1,_mm_setzero_si128()),ly1);
        	C0= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp1,_mm_setzero_si128()),ly1);

        	c1= _mm_mullo_epi16(_mm_unpacklo_epi8(tmp2,_mm_setzero_si128()),t0);
        	C1= _mm_mullo_epi16(_mm_unpackhi_epi8(tmp2,_mm_setzero_si128()),t0);

        	c2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpacklo_epi8(tx,_mm_setzero_si128()));
        	C2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpackhi_epi8(tx,_mm_setzero_si128()));



        	c0= _mm_add_epi16(c0,c1);
        	c2= _mm_add_epi16(c2,c3);
        	c2= _mm_add_epi16(c2,add);
        	c0= _mm_add_epi16(c0,c2);

        	C0= _mm_add_epi16(C0,C1);
        	C2= _mm_add_epi16(C2,c3);
        	C2= _mm_add_epi16(C2,add);
        	C0= _mm_add_epi16(C0,C2);

        	c0= _mm_srli_epi16(c0,6);
        	C0= _mm_srli_epi16(C0,6);

        	c0= _mm_packus_epi16(c0,C0);

        	_mm_storeu_si128((__m128i*)(src + y*stride), c0);

        	// second half of 32

        	c0= _mm_mullo_epi16(_mm_unpacklo_epi8(TMP1,_mm_setzero_si128()),ly1);
        	C0= _mm_mullo_epi16(_mm_unpackhi_epi8(TMP1,_mm_setzero_si128()),ly1);

        	c1= _mm_mullo_epi16(_mm_unpacklo_epi8(TMP2,_mm_setzero_si128()),t0);
        	C1= _mm_mullo_epi16(_mm_unpackhi_epi8(TMP2,_mm_setzero_si128()),t0);

        	c2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpacklo_epi8(TX,_mm_setzero_si128()));
        	C2= _mm_mullo_epi16(_mm_set1_epi16(31 - y),_mm_unpackhi_epi8(TX,_mm_setzero_si128()));



        	c0= _mm_add_epi16(c0,c1);
        	c2= _mm_add_epi16(c2,c3);
        	c2= _mm_add_epi16(c2,add);
        	c0= _mm_add_epi16(c0,c2);

        	C0= _mm_add_epi16(C0,C1);
        	C2= _mm_add_epi16(C2,c3);
        	C2= _mm_add_epi16(C2,add);
        	C0= _mm_add_epi16(C0,C2);

        	c0= _mm_srli_epi16(c0,6);
        	C0= _mm_srli_epi16(C0,6);

        	c0= _mm_packus_epi16(c0,C0);

        	_mm_storeu_si128((__m128i*)(src + 16 + y*stride), c0);

        	LY= _mm_srli_si128(LY,1);
        }


}
#else
static void FUNCC(pred_planar_0)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    for (y = 0; y < 4; y++)
        for (x = 0; x < 4; x++)
            POS(x, y) = ((3 - x) * left[y]  + (x + 1) * top[4] +
                         (3 - y) * top[x] + (y + 1) * left[4] + 4) >>
                        (3);
}
static void FUNCC(pred_planar_1)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    for (y = 0; y < 8; y++)
        for (x = 0; x < 8; x++)
            POS(x, y) = ((7 - x) * left[y]  + (x + 1) * top[8] +
                         (7 - y) * top[x] + (y + 1) * left[8] + 8) >>
                        (4);
}
static void FUNCC(pred_planar_2)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    for (y = 0; y < 16; y++)
        for (x = 0; x < 16; x++)
            POS(x, y) = ((15 - x) * left[y]  + (x + 1) * top[16] +
                         (15 - y) * top[x] + (y + 1) * left[16] + 16) >>
                        (5);
}
static void FUNCC(pred_planar_3)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                               ptrdiff_t stride, int log2_size)
{
    int x, y;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    for (y = 0; y < 32; y++)
        for (x = 0; x < 32; x++)
            POS(x, y) = ((31 - x) * left[y]  + (x + 1) * top[32] +
                         (31 - y) * top[x] + (y + 1) * left[32] + 32) >>
                        (6);
}
#endif

static void FUNCC(pred_dc)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                           ptrdiff_t stride, int log2_size, int c_idx)
{
    int i, j, x, y;
    int size = (1 << log2_size);
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;
    int dc = size;
    pixel4 a;
    for (i = 0; i < size; i++)
        dc += left[i] + top[i];

    dc >>= log2_size + 1;

    a = PIXEL_SPLAT_X4(dc);

    for (i = 0; i < size; i++)
        for (j = 0; j < size / 4; j++)
            AV_WN4PA(&POS(j * 4, i), a);

    if (c_idx == 0 && size < 32) {
        POS(0, 0) = (left[0] + 2 * dc  + top[0] + 2) >> 2;
        for (x = 1; x < size; x++)
            POS(x, 0) = (top[x] + 3 * dc + 2) >> 2;
        for (y = 1; y < size; y++)
            POS(0, y) = (left[y] + 3 * dc + 2) >> 2;
    }
}

static void FUNCC(pred_angular)(uint8_t *_src, const uint8_t *_top, const uint8_t *_left,
                                ptrdiff_t stride, int log2_size, int c_idx, int mode)
{
    int x, y;
    int size = 1 << log2_size;
    pixel *src = (pixel*)_src;
    const pixel *top = (const pixel*)_top;
    const pixel *left = (const pixel*)_left;

    const int intra_pred_angle[] = {
        32, 26, 21, 17, 13, 9, 5, 2, 0, -2, -5, -9, -13, -17, -21, -26, -32,
        -26, -21, -17, -13, -9, -5, -2, 0, 2, 5, 9, 13, 17, 21, 26, 32
    };
    const int inv_angle[] = {
        -4096, -1638, -910, -630, -482, -390, -315, -256, -315, -390, -482,
        -630, -910, -1638, -4096
    };

    int angle = intra_pred_angle[mode-2];
    pixel ref_array[3*MAX_TB_SIZE+1];
    const pixel *ref;
    int last = (size * angle) >> 5;

    if (mode >= 18) {
        ref = top - 1;
        if (angle < 0 && last < -1) {
            for (x = 0; x <= size; x++)
                (ref_array + size)[x] = top[x - 1];
            for (x = last; x <= -1; x++)
                (ref_array + size)[x] = left[-1 + ((x * inv_angle[mode-11] + 128) >> 8)];
            ref = ref_array + size;
        }

        for (y = 0; y < size; y++) {
            int idx = ((y + 1) * angle) >> 5;
            int fact = ((y + 1) * angle) & 31;
            if (fact) {
                for (x = 0; x < size; x++) {
                    POS(x, y) = ((32 - fact) * ref[x + idx + 1] + fact * ref[x + idx + 2] + 16) >> 5;
                }
            } else {
                for (x = 0; x < size; x++) {
                    POS(x, y) = ref[x + idx + 1];
                }
            }
        }
        if (mode == 26 && c_idx == 0 && size < 32) {
            for (y = 0; y < size; y++)
                POS(0, y) = av_clip_pixel(top[0] + ((left[y] - left[-1]) >> 1));
        }
    } else {
        ref = left - 1;
        if (angle < 0 && last < -1) {
            for (x = 0; x <= size; x++)
                (ref_array + size)[x] = left[x - 1];
            for (x = last; x <= -1; x++)
                (ref_array + size)[x] = top[-1 + ((x * inv_angle[mode-11] + 128) >> 8)];
            ref = ref_array + size;
        }

        for (x = 0; x < size; x++) {
            int idx = ((x + 1) * angle) >> 5;
            int fact = ((x + 1) * angle) & 31;
            if (fact) {
                for (y = 0; y < size; y++) {
                    POS(x, y) = ((32 - fact) * ref[y + idx + 1] + fact * ref[y + idx + 2] + 16) >> 5;
                }
            } else {
                for (y = 0; y < size; y++) {
                    POS(x, y) = ref[y + idx + 1];
                }
            }
        }
        if (mode == 10 && c_idx == 0 && size < 32) {
            for (x = 0; x < size; x++)
                POS(x, 0) = av_clip_pixel(left[0] + ((top[x] - top[-1]) >> 1));
        }
    }
}

#undef POS
