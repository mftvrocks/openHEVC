/*
 * Copyright (c) 2013 Seppo Tomperi
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

#include "config.h"
#include "libavutil/cpu.h"
#include "libavutil/x86/asm.h"
#include "libavutil/x86/cpu.h"
#include "libavcodec/get_bits.h" /* required for hevcdsp.h GetBitContext */
#include "libavcodec/hevcdsp.h"


/***********************************/
/* deblocking */

#define LFC_FUNC(DIR, DEPTH, OPT)                                        \
void ff_hevc_ ## DIR ## _loop_filter_chroma_ ## DEPTH ## _ ## OPT(uint8_t *_pix, ptrdiff_t _stride, int *_tc, uint8_t *_no_p, uint8_t *_no_q);

#define LFL_FUNC(DIR, DEPTH, OPT)                                        \
void ff_hevc_ ## DIR ## _loop_filter_luma_ ## DEPTH ## _ ## OPT(uint8_t *_pix, ptrdiff_t stride, int *_beta, int *_tc, \
uint8_t *_no_p, uint8_t *_no_q);

#define LFC_FUNCS(type, depth) \
LFC_FUNC(h, depth, sse2)  \
LFC_FUNC(v, depth, sse2)

#define LFL_FUNCS(type, depth) \
LFL_FUNC(h, depth, ssse3)  \
LFL_FUNC(v, depth, ssse3)


LFC_FUNCS(uint8_t,   8)
LFL_FUNCS(uint8_t,   8)

void ff_hevc_dequant4x4_sse4(int16_t *coeffs, int qp);
void ff_hevc_dequant8x8_sse4(int16_t *coeffs, int qp);
void ff_hevc_dequant16x16_sse4(int16_t *coeffs, int qp);
void ff_hevc_dequant32x32_sse4(int16_t *coeffs, int qp);

//IDCT functions
void ff_hevc_transform_4x4_luma_add_8_sse4(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);

void ff_hevc_transform_4x4_add_8_sse4(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);
void ff_hevc_transform_8x8_add_8_sse4(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);
void ff_hevc_transform_16x16_add_8_sse4(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);
void ff_hevc_transform_32x32_add_8_sse4(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride);

// MC functions
void ff_hevc_put_unweighted_pred_sse(uint8_t *_dst, ptrdiff_t _dststride,int16_t *src, ptrdiff_t srcstride,int width, int height);
void ff_hevc_put_weighted_pred_avg_sse(uint8_t *_dst, ptrdiff_t _dststride,int16_t *src1, int16_t *src2, ptrdiff_t srcstride,int width, int height);
void ff_hevc_weighted_pred_sse(uint8_t denom, int16_t wlxFlag, int16_t olxFlag,uint8_t *_dst, ptrdiff_t _dststride,int16_t *src, ptrdiff_t srcstride,int width, int height);
void ff_hevc_weighted_pred_avg_sse(uint8_t denom, int16_t wl0Flag, int16_t wl1Flag,int16_t ol0Flag, int16_t ol1Flag, uint8_t *_dst, ptrdiff_t _dststride,int16_t *src1, int16_t *src2, ptrdiff_t srcstride,int width, int height);

void ff_hevc_put_hevc_epel_pixels_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);
void ff_hevc_put_hevc_epel_pixels_8_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);

void ff_hevc_put_hevc_epel_h_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);
void ff_hevc_put_hevc_epel_h_8_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);

void ff_hevc_put_hevc_epel_v_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);
void ff_hevc_put_hevc_epel_v_8_sse(int16_t *dst, ptrdiff_t dststride, uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);

void ff_hevc_put_hevc_epel_hv_8_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int mx, int my, int16_t* mcbuffer);

void ff_hevc_put_hevc_qpel_pixels_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_1_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_2_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_3_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_v_1_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_v_2_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_v_3_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);

void ff_hevc_put_hevc_qpel_h_1_v_1_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_1_v_2_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_1_v_3_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_2_v_1_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_2_v_2_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_2_v_3_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_3_v_1_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_3_v_2_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);
void ff_hevc_put_hevc_qpel_h_3_v_3_sse(int16_t *dst, ptrdiff_t dststride,uint8_t *_src, ptrdiff_t _srcstride,int width, int height, int16_t* mcbuffer);

// SAO functions

void ff_hevc_sao_edge_filter_0_8_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx);
void ff_hevc_sao_edge_filter_1_8_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx);
void ff_hevc_sao_edge_filter_2_8_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx);
void ff_hevc_sao_edge_filter_3_8_sse(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx);


//LF_FUNCS(uint16_t, 10)


void ff_hevcdsp_init_x86(HEVCDSPContext *c, const int bit_depth, const int pcm_deblock)
{
    int mm_flags = av_get_cpu_flags();

    if (bit_depth == 8) {
        if (EXTERNAL_MMX(mm_flags)) {
            /*if (mm_flags & AV_CPU_FLAG_CMOV)
                c->h264_luma_dc_dequant_idct = ff_h264_luma_dc_dequant_idct_mmx; */

            if (EXTERNAL_MMXEXT(mm_flags)) {
#if ARCH_X86_32 && HAVE_MMXEXT_EXTERNAL
                /* MMEXT optimizations */
#endif /* ARCH_X86_32 && HAVE_MMXEXT_EXTERNAL */

                if (EXTERNAL_SSE2(mm_flags)) {
                    if (!pcm_deblock) {
                        c->hevc_v_loop_filter_chroma = ff_hevc_v_loop_filter_chroma_8_sse2;
                        c->hevc_h_loop_filter_chroma = ff_hevc_h_loop_filter_chroma_8_sse2;
                    }
                }
                if (EXTERNAL_SSSE3(mm_flags)) {
#if ARCH_X86_64
                        c->hevc_v_loop_filter_luma = ff_hevc_v_loop_filter_luma_8_ssse3;
                        c->hevc_h_loop_filter_luma = ff_hevc_h_loop_filter_luma_8_ssse3;
#endif

                }
                if (EXTERNAL_SSE4(mm_flags)) {
#if ARCH_X86_64
#if GCC_VERSION > MIN_GCC_VERSION_MC
                	c->put_unweighted_pred = ff_hevc_put_unweighted_pred_sse;
                	c->put_weighted_pred_avg = ff_hevc_put_weighted_pred_avg_sse;
                	c->weighted_pred = ff_hevc_weighted_pred_sse;
                	c->weighted_pred_avg = ff_hevc_weighted_pred_avg_sse;


                	c->put_hevc_epel[0][0] = ff_hevc_put_hevc_epel_pixels_8_sse;
                	c->put_hevc_epel[0][1] = ff_hevc_put_hevc_epel_h_8_sse;
                	c->put_hevc_epel[1][0] = ff_hevc_put_hevc_epel_v_8_sse;
                	c->put_hevc_epel[1][1] = ff_hevc_put_hevc_epel_hv_8_sse;

                	c->put_hevc_qpel[0][0] = ff_hevc_put_hevc_qpel_pixels_sse;
                	c->put_hevc_qpel[0][1] = ff_hevc_put_hevc_qpel_h_1_sse;
                	c->put_hevc_qpel[0][2] = ff_hevc_put_hevc_qpel_h_2_sse;
                	c->put_hevc_qpel[0][3] = ff_hevc_put_hevc_qpel_h_3_sse;
                	c->put_hevc_qpel[1][0] = ff_hevc_put_hevc_qpel_v_1_sse;
                	c->put_hevc_qpel[2][0] = ff_hevc_put_hevc_qpel_v_2_sse;
                	c->put_hevc_qpel[3][0] = ff_hevc_put_hevc_qpel_v_3_sse;

                	c->put_hevc_qpel[1][1] = ff_hevc_put_hevc_qpel_h_1_v_1_sse;
                	c->put_hevc_qpel[1][2] = ff_hevc_put_hevc_qpel_h_2_v_1_sse;
                	c->put_hevc_qpel[1][3] = ff_hevc_put_hevc_qpel_h_3_v_1_sse;
                	c->put_hevc_qpel[2][1] = ff_hevc_put_hevc_qpel_h_1_v_2_sse;
                	c->put_hevc_qpel[2][2] = ff_hevc_put_hevc_qpel_h_2_v_2_sse;
                	c->put_hevc_qpel[2][3] = ff_hevc_put_hevc_qpel_h_3_v_2_sse;
                	c->put_hevc_qpel[3][1] = ff_hevc_put_hevc_qpel_h_1_v_3_sse;
                	c->put_hevc_qpel[3][2] = ff_hevc_put_hevc_qpel_h_2_v_3_sse;
                	c->put_hevc_qpel[3][3] = ff_hevc_put_hevc_qpel_h_3_v_3_sse;
#endif
#endif


                	c->sao_edge_filter[0] = ff_hevc_sao_edge_filter_0_8_sse;
                	c->sao_edge_filter[1] = ff_hevc_sao_edge_filter_1_8_sse;
                	c->sao_edge_filter[2] = ff_hevc_sao_edge_filter_2_8_sse;
                	c->sao_edge_filter[3] = ff_hevc_sao_edge_filter_3_8_sse;


                    c->dequant[0] = ff_hevc_dequant4x4_sse4;
                    c->dequant[1] = ff_hevc_dequant8x8_sse4;
                    c->dequant[2] = ff_hevc_dequant16x16_sse4;
                    c->dequant[3] = ff_hevc_dequant32x32_sse4;

                    c->transform_4x4_luma_add = ff_hevc_transform_4x4_luma_add_8_sse4;

                    c->transform_add[0] = ff_hevc_transform_4x4_add_8_sse4;
                    c->transform_add[1] = ff_hevc_transform_8x8_add_8_sse4;
                    c->transform_add[2] = ff_hevc_transform_16x16_add_8_sse4;
                    c->transform_add[3] = ff_hevc_transform_32x32_add_8_sse4;
                    }
                if (EXTERNAL_AVX(mm_flags)) {

                }
            }
        }
    } else if (bit_depth == 10) {
        if (EXTERNAL_MMX(mm_flags)) {
            if (EXTERNAL_MMXEXT(mm_flags)) {
#if ARCH_X86_32
#endif /* ARCH_X86_32 */
                if (EXTERNAL_SSE2(mm_flags)) {
#if HAVE_ALIGNED_STACK
                    /*stuff that requires aligned stack */
#endif /* HAVE_ALIGNED_STACK */
                }
                if (EXTERNAL_SSE4(mm_flags)) {
                }
                if (EXTERNAL_AVX(mm_flags)) {
                }
            }
        }
    }
}
