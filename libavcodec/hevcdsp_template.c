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

#include "libavutil/avassert.h"
#include "libavutil/pixdesc.h"
#include "get_bits.h"
#include "bit_depth_template.c"
#include "hevcdata.h"
#include "hevcdsp.h"
#include "hevc.h"
#define shift_1st 7
#define add_1st (1 << (shift_1st - 1))
#define shift_2nd (20 - BIT_DEPTH)
#define add_2nd (1 << (shift_2nd - 1))

#if __GNUC__
//#define GCC_OPTIMIZATION_ENABLE
#endif
#define OPTIMIZATION_ENABLE



#define SET(dst, x) (dst) = (x)
#define SCALE(dst, x) (dst) = av_clip_int16(((x) + add) >> shift)
#define ADD_AND_SCALE(dst, x) (dst) = av_clip_pixel((dst) + av_clip_int16(((x) + add) >> shift))

static void FUNC(put_pcm)(uint8_t *_dst, ptrdiff_t _stride, int size,
                          GetBitContext *gb, int pcm_bit_depth)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);

    for (y = 0; y < size; y++) {
        for (x = 0; x < size; x++)
            dst[x] = get_bits(gb, pcm_bit_depth) << (BIT_DEPTH - pcm_bit_depth);
        dst += stride;
    }
}

static void FUNC(dequant4x4)(int16_t *coeffs, int qp)
{
    int x, y;
    int size = 4;

    const uint8_t level_scale[] = { 40, 45, 51, 57, 64, 72 };

    //TODO: scaling_list_enabled_flag support

    int shift  = BIT_DEPTH -3;
    int scale  = level_scale[qp % 6] << (qp / 6);
    int add    = 1 << (shift - 1);
    int scale2 = scale << 4;
            for (y = 0; y < 4*4; y+=4) {
                SCALE(coeffs[y+0], (coeffs[y+0] * scale2));
                SCALE(coeffs[y+1], (coeffs[y+1] * scale2));
                SCALE(coeffs[y+2], (coeffs[y+2] * scale2));
                SCALE(coeffs[y+3], (coeffs[y+3] * scale2));
            }

}
static void FUNC(dequant8x8)(int16_t *coeffs, int qp)
{
    int x, y;
    int size = 8;

    const uint8_t level_scale[] = { 40, 45, 51, 57, 64, 72 };

    //TODO: scaling_list_enabled_flag support

    int shift  = BIT_DEPTH - 2;
    int scale  = level_scale[qp % 6] << (qp / 6);
    int add    = 1 << (shift - 1);
    int scale2 = scale << 4;
#define FILTER()                                    \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2))

            for (y = 0; y < 8*8; y+=8) {
                for (x = 0; x < 8; x++) {
                    FILTER();
                }
            }
#undef FILTER
}

static void FUNC(dequant16x16)(int16_t *coeffs, int qp)
{
    int x, y;
    int size = 16;

    const uint8_t level_scale[] = { 40, 45, 51, 57, 64, 72 };

    //TODO: scaling_list_enabled_flag support

    int shift  = BIT_DEPTH -1;
    int scale  = level_scale[qp % 6] << (qp / 6);
    int add    = 1 << (shift - 1);
    int scale2 = scale << 4;
#define FILTER()                                    \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2))
            for (y = 0; y < 16*16; y+=16) {
                for (x = 0; x < 16; x++) {
                    FILTER();
                }
            }
#undef FILTER
}

static void FUNC(dequant32x32)(int16_t *coeffs, int qp)
{
    int x, y;
    int size = 32;

    const uint8_t level_scale[] = { 40, 45, 51, 57, 64, 72 };

    //TODO: scaling_list_enabled_flag support

    int shift  = BIT_DEPTH;
    int scale  = level_scale[qp % 6] << (qp / 6);
    int add    = 1 << (shift - 1);
    int scale2 = scale << 4;
#define FILTER()                                    \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2)); \
        x++;                                        \
        SCALE(coeffs[y+x], (coeffs[y+x] * scale2))
            for (y = 0; y < 32*32; y+=32) {
                for (x = 0; x < 32; x++) {
                    FILTER();
                }
            }

#undef FILTER
}
static void FUNC(transquant_bypass4x4)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int size = 1 << 2;

            for (y = 0; y < 4; y++) {
                for (x = 0; x < 4; x++) {
                    dst[x] += *coeffs;
                    coeffs++;
                }
                dst += stride;
            }

}

static void FUNC(transquant_bypass8x8)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int size = 1 << 3;

            for (y = 0; y < 8; y++) {
                for (x = 0; x < 8; x++) {
                    dst[x] += *coeffs;
                    coeffs++;
                }
                dst += stride;
            }
}

static void FUNC(transquant_bypass16x16)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int size = 1 << 4;

            for (y = 0; y < 16; y++) {
                for (x = 0; x < 16; x++) {
                    dst[x] += *coeffs;
                    coeffs++;
                }
                dst += stride;
            }

}

static void FUNC(transquant_bypass32x32)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int size = 1 << 5;

            for (y = 0; y < 32; y++) {
                for (x = 0; x < 32; x++) {
                    dst[x] += *coeffs;
                    coeffs++;
                }
                dst += stride;
            }

}



static void FUNC(transform_skip)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int size = 4;
    int shift = 13 - BIT_DEPTH;
#if BIT_DEPTH <= 13
    int offset = 1 << (shift - 1);
#else
    int offset = 0;
#endif
    int x, y;
    switch (size){
        case 32:
            for (y = 0; y < 32*32; y+=32) {
                for (x = 0; x < 32; x++) {
                    dst[x] = av_clip_pixel(dst[x] + ((coeffs[y + x] + offset) >> shift));
                }
                dst += stride;
            }
            break;
        case 16:
            for (y = 0; y < 16*16; y+=16) {
                for (x = 0; x < 16; x++) {
                    dst[x] = av_clip_pixel(dst[x] + ((coeffs[y + x] + offset) >> shift));
                }
                dst += stride;
            }
            break;
        case 8:
            for (y = 0; y < 8*8; y+=8) {
                for (x = 0; x < 8; x++) {
                    dst[x] = av_clip_pixel(dst[x] + ((coeffs[y + x] + offset) >> shift));
                }
                dst += stride;
            }
            break;
        case 4:
            for (y = 0; y < 4*4; y+=4) {
                for (x = 0; x < 4; x++) {
                    dst[x] = av_clip_pixel(dst[x] + ((coeffs[y + x] + offset) >> shift));
                }
                dst += stride;
            }
            break;
    }
#undef FILTER
}

static void FUNC(transform_4x4_luma_add)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int i;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int16_t *src = coeffs;

    
#define TR_4x4_LUMA(dst, src, step, assign)                                     \
do {                                                                        \
int c0 = src[0*step] + src[2*step];                                     \
int c1 = src[2*step] + src[3*step];                                     \
int c2 = src[0*step] - src[3*step];                                     \
int c3 = 74 * src[1*step];                                              \
\
assign(dst[2*step], 74 * (src[0*step] - src[2*step] + src[3*step]));    \
assign(dst[0*step], 29 * c0 + 55 * c1 + c3);                            \
assign(dst[1*step], 55 * c2 - 29 * c1 + c3);                            \
assign(dst[3*step], 55 * c0 + 29 * c2 - c3);                            \
} while (0)
    
    int shift = 7;
    int add = 1 << (shift - 1);
    for (i = 0; i < 4; i++) {
        TR_4x4_LUMA(src, src, 4, SCALE);
        src++;
    }
    
    shift = 20 - BIT_DEPTH;
    add = 1 << (shift - 1);
    for (i = 0; i < 4; i++) {
        TR_4x4_LUMA(dst, coeffs, 1, ADD_AND_SCALE);
        coeffs += 4;
        dst += stride;
    }
    
#undef TR_4x4_LUMA

}

#define TR_4(dst, src, dstep, sstep, assign)                                    \
    do {                                                                        \
        const int e0 = transform[8*0][0] * src[0*sstep] +                       \
                       transform[8*2][0] * src[2*sstep];                        \
        const int e1 = transform[8*0][1] * src[0*sstep] +                       \
                       transform[8*2][1] * src[2*sstep];                        \
        const int o0 = transform[8*1][0] * src[1*sstep] +                       \
                       transform[8*3][0] * src[3*sstep];                        \
        const int o1 = transform[8*1][1] * src[1*sstep] +                       \
                       transform[8*3][1] * src[3*sstep];                        \
                                                                                \
        assign(dst[0*dstep], e0 + o0);                                          \
        assign(dst[1*dstep], e1 + o1);                                          \
        assign(dst[2*dstep], e1 - o1);                                          \
        assign(dst[3*dstep], e0 - o0);                                          \
    } while (0)
#define TR_4_1(dst, src) TR_4(dst, src, 4, 4, SCALE)
#define TR_4_2(dst, src) TR_4(dst, src, 1, 1, ADD_AND_SCALE)

static void FUNC(transform_4x4_add)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int i;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int16_t *src = coeffs;

    int shift = 7;
    int add = 1 << (shift - 1);
    for (i = 0; i < 4; i++) {
        TR_4_1(src, src);
        src++;
    }
    
    shift = 20 - BIT_DEPTH;
    add = 1 << (shift - 1);
    for (i = 0; i < 4; i++) {
        TR_4_2(dst, coeffs);
        coeffs += 4;
        dst += stride;
    }

}

#define TR_8(dst, src, dstep, sstep, assign)                \
    do {                                                    \
        int i, j;                                           \
        int e_8[4];                                         \
        int o_8[4] = { 0 };                                 \
        for (i = 0; i < 4; i++)                             \
            for (j = 1; j < 8; j += 2)                      \
                o_8[i] += transform[4*j][i] * src[j*sstep]; \
        TR_4(e_8, src, 1, 2*sstep, SET);                    \
for (i = 0; i < 4; i++) {                           \
            assign(dst[i*dstep], e_8[i] + o_8[i]);          \
            assign(dst[(7-i)*dstep], e_8[i] - o_8[i]);      \
        }                                                   \
    } while (0)
#define TR_16(dst, src, dstep, sstep, assign)                   \
    do {                                                        \
        int i, j;                                               \
        int e_16[8];                                            \
        int o_16[8] = { 0 };                                    \
        for (i = 0; i < 8; i++)                                 \
            for (j = 1; j < 16; j += 2)                         \
                o_16[i] += transform[2*j][i] * src[j*sstep];    \
        TR_8(e_16, src, 1, 2*sstep, SET);                       \
                                                                 \
        for (i = 0; i < 8; i++) {                               \
            assign(dst[i*dstep], e_16[i] + o_16[i]);            \
            assign(dst[(15-i)*dstep], e_16[i] - o_16[i]);       \
        }                                                       \
    } while (0)
#define TR_32(dst, src, dstep, sstep, assign)               \
    do {                                                    \
        int i, j;                                           \
        int e_32[16];                                       \
        int o_32[16] = { 0 };                               \
        for (i = 0; i < 16; i++)                            \
            for (j = 1; j < 32; j += 2)                     \
                o_32[i] += transform[j][i] * src[j*sstep];  \
        TR_16(e_32, src, 1, 2*sstep, SET);                \
                                                            \
        for (i = 0; i < 16; i++) {                          \
            assign(dst[i*dstep], e_32[i] + o_32[i]);        \
            assign(dst[(31-i)*dstep], e_32[i] - o_32[i]);   \
        }                                                   \
    } while (0)

#define TR_8_1(dst, src) TR_8(dst, src, 8, 8, SCALE)
#define TR_16_1(dst, src) TR_16(dst, src, 16, 16, SCALE)
#define TR_32_1(dst, src) TR_32(dst, src, 32, 32, SCALE)

#define TR_8_2(dst, src) TR_8(dst, src, 1, 1, ADD_AND_SCALE)
#define TR_16_2(dst, src) TR_16(dst, src, 1, 1, ADD_AND_SCALE)
#define TR_32_2(dst, src) TR_32(dst, src, 1, 1, ADD_AND_SCALE)

static void FUNC(transform_8x8_add)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride) {
    int i;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int16_t *src = coeffs;

    int shift = 7;
    int add = 1 << (shift - 1);
    for (i = 0; i < 8; i++) {
        TR_8_1(src, src);
        src++;
    }
    src = coeffs;
    shift = 20 - BIT_DEPTH;
    add = 1 << (shift - 1);
    for (i = 0; i < 8; i++) {
        TR_8_2(dst, coeffs);
        coeffs += 8;
        dst += stride;
    }

}

static void FUNC(transform_16x16_add)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{
    int i;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int16_t *src = coeffs;
    int32_t shift;

    shift = 7;
    int add = 1 << (shift - 1);
    for (i = 0; i < 16; i++) {
        TR_16_1(src, src);
        src++;
    }
    
    shift = 20 - BIT_DEPTH;
    add = 1 << (shift - 1);

    for (i = 0; i < 16; i++) {
        TR_16_2(dst, coeffs);
        coeffs += 16;
        dst += stride;
    }
}

static void FUNC(transform_32x32_add)(uint8_t *_dst, int16_t *coeffs, ptrdiff_t _stride)
{

#define IT32x32_even(i,w) ( src[ 0*w] * transform[ 0][i] ) + ( src[16*w] * transform[16][i] )
#define IT32x32_odd(i,w)  ( src[ 8*w] * transform[ 8][i] ) + ( src[24*w] * transform[24][i] )
#define IT16x16(i,w)      ( src[ 4*w] * transform[ 4][i] ) + ( src[12*w] * transform[12][i] ) + ( src[20*w] * transform[20][i] ) + ( src[28*w] * transform[28][i] )
#define IT8x8(i,w)        ( src[ 2*w] * transform[ 2][i] ) + ( src[ 6*w] * transform[ 6][i] ) + ( src[10*w] * transform[10][i] ) + ( src[14*w] * transform[14][i] ) + \
( src[18*w] * transform[18][i] ) + ( src[22*w] * transform[22][i] ) + ( src[26*w] * transform[26][i] ) + ( src[30*w] * transform[30][i] )
#define IT4x4(i,w)        ( src[ 1*w] * transform[ 1][i] ) + ( src[ 3*w] * transform[ 3][i] ) + ( src[ 5*w] * transform[ 5][i] ) + ( src[ 7*w] * transform[ 7][i] ) + \
( src[ 9*w] * transform[ 9][i] ) + ( src[11*w] * transform[11][i] ) + ( src[13*w] * transform[13][i] ) + ( src[15*w] * transform[15][i] ) + \
( src[17*w] * transform[17][i] ) + ( src[19*w] * transform[19][i] ) + ( src[21*w] * transform[21][i] ) + ( src[23*w] * transform[23][i] ) + \
( src[25*w] * transform[25][i] ) + ( src[27*w] * transform[27][i] ) + ( src[29*w] * transform[29][i] ) + ( src[31*w] * transform[31][i] )
    int i,j;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t stride = _stride / sizeof(pixel);
    int shift = 7;
    int add = 1 << (shift - 1);
    int16_t *src = coeffs;
    for (i = 0; i < 32; i++) {
#ifndef OPTIMIZATION_ENABLE
        TR_32_1(src, src);
#else
        int e_32[16];
        int o_32[16];
        int e_16[8];
        int e_8[4];
        int odd;
        const int e0 = IT32x32_even(0,32);
        const int e1 = IT32x32_even(1,32);
        const int o0 = IT32x32_odd(0,32);
        const int o1 = IT32x32_odd(1,32);
        e_8[0] = e0 + o0;
        e_8[1] = e1 + o1;
        e_8[2] = e1 - o1;
        e_8[3] = e0 - o0;
        for (j = 0; j < 4; j++) {
            odd       = IT16x16(j,32);
            e_16[  j] = e_8[j] + odd;
            e_16[7-j] = e_8[j] - odd;
        }
        for (j = 0; j < 8; j++) {
            odd        = IT8x8(j,32);
            e_32[   j] = e_16[j] + odd;
            e_32[15-j] = e_16[j] - odd;
        }
        for (j = 0; j < 16; j++){
            o_32[j] = IT4x4(j,32);
            
        }
        
        for (j = 0; j < 16; j++) {
            odd        = o_32[j];
            
            SCALE(src[(   j)*32], (e_32[j] + odd));
            SCALE(src[(31-j)*32], (e_32[j] - odd));
        }
        
#endif
        src++;
    }
    src   = coeffs;
    shift = 20 - BIT_DEPTH;
    add   = 1 << (shift - 1);
    for (i = 0; i < 32; i++) {
#ifndef OPTIMIZATION_ENABLE
        TR_32_2(dst, coeffs);
        coeffs += 32;
#else
        int e_32[16];
        int e_16[8];
        int e_8[4];
        int odd;
        const int e0 = IT32x32_even(0,1);
        const int e1 = IT32x32_even(1,1);
        const int o0 = IT32x32_odd(0,1);
        const int o1 = IT32x32_odd(1,1);
        e_8[0] = e0 + o0;
        e_8[1] = e1 + o1;
        e_8[2] = e1 - o1;
        e_8[3] = e0 - o0;
        for (j = 0; j < 4; j++) {
            odd       = IT16x16(j,1);
            e_16[  j] = e_8[j] + odd;
            e_16[7-j] = e_8[j] - odd;
        }
        for (j = 0; j < 8; j++) {
            odd        = IT8x8(j,1);
            e_32[   j] = e_16[j] + odd;
            e_32[15-j] = e_16[j] - odd;
        }
        pixel dst1[32*32];
        memcpy(dst1, dst, 32*32*sizeof(pixel));
        for (j = 0; j < 16; j++) {
            odd       = IT4x4(j,1);
            ADD_AND_SCALE(dst[   j], e_32[j] + odd);
            ADD_AND_SCALE(dst[31-j], e_32[j] - odd);
        }
        src += 32;
#endif
        dst += stride;
    }
#undef IT32x32_even
#undef IT32x32_odd
#undef IT16x16
#undef IT8x8
#undef IT4x4
}


static void FUNC(sao_band_filter_wpp_0)( uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int width, int height, int c_idx)
{
    uint8_t *dst = _dst;
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int band_table[32] = { 0 };
    int k, y, x;
    int chroma = c_idx!=0;
    int shift = BIT_DEPTH - 5;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_left_class = sao->band_position[c_idx];

    int init_y = 0, init_x =0;

            if(!borders[2] )
                width -= ((8>>chroma)+2) ;
            if(!borders[3] )
                height -= ((4>>chroma)+2);

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    for (k = 0; k < 4; k++)
        band_table[(k + sao_left_class) & 31] = k + 1;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
                x++;
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
            }
            dst += stride;
            src += stride;
        }
}

static void FUNC(sao_band_filter_wpp_1)( uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int width, int height, int c_idx)
{
    uint8_t *dst = _dst;
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int band_table[32] = { 0 };
    int k, y, x;
    int chroma = c_idx!=0;
    int shift = BIT_DEPTH - 5;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_left_class = sao->band_position[c_idx];

    int init_y = 0, init_x =0;


            init_y = -(4>>chroma)-2;
            if(!borders[2] )
                width -= ((8>>chroma)+2);
            height = (4>>chroma)+2;

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    for (k = 0; k < 4; k++)
        band_table[(k + sao_left_class) & 31] = k + 1;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
                x++;
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
            }
            dst += stride;
            src += stride;
        }
}

static void FUNC(sao_band_filter_wpp_2)( uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int width, int height, int c_idx)
{
    uint8_t *dst = _dst;
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int band_table[32] = { 0 };
    int k, y, x;
    int chroma = c_idx!=0;
    int shift = BIT_DEPTH - 5;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_left_class = sao->band_position[c_idx];

    int init_y = 0, init_x =0;


            init_x = -(8>>chroma)-2;
            width = (8>>chroma)+2;
            if(!borders[3])
                height -= ((4>>chroma)+2);

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    for (k = 0; k < 4; k++)
        band_table[(k + sao_left_class) & 31] = k + 1;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
                x++;
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
            }
            dst += stride;
            src += stride;
        }
}

static void FUNC(sao_band_filter_wpp_3)( uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int width, int height, int c_idx)
{
    uint8_t *dst = _dst;
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int band_table[32] = { 0 };
    int k, y, x;
    int chroma = c_idx!=0;
    int shift = BIT_DEPTH - 5;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_left_class = sao->band_position[c_idx];

    int init_y = 0, init_x =0;

            init_y = -(4>>chroma)-2;
            init_x = -(8>>chroma)-2;
            width = (8>>chroma)+2;
            height = (4>>chroma)+2;

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    for (k = 0; k < 4; k++)
        band_table[(k + sao_left_class) & 31] = k + 1;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
                x++;
                dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
            }
            dst += stride;
            src += stride;
        }
}


static void FUNC(sao_edge_filter_wpp_0)(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx)
{
    int x, y;
    uint8_t *dst = _dst;   // put here pixel
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int chroma = c_idx!=0;
    //struct SAOParams *sao;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_eo_class = sao->eo_class[c_idx];

    const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };
    const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };

    int init_x = 0, init_y = 0, width = _width, height = _height;

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))

            if(!borders[2] )
                width -= ((8>>chroma)+2) ;
            if(!borders[3] )
                height -= ((4>>chroma)+2);

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    init_y = init_x = 0;
    if (sao_eo_class != SAO_EO_VERT) {
        if (borders[0]) {
            int offset_val = sao_offset_val[0];
            int y_stride   = 0;
            for (y = 0; y < height; y++) {
                dst[y_stride] = av_clip_pixel(src[y_stride] + offset_val);
                y_stride += stride;
            }
            init_x = 1;
        }
        if (borders[2]) {
            int offset_val = sao_offset_val[0];
            int x_stride   = _width-1;
            for (x = 0; x < height; x++) {
                dst[x_stride] = av_clip_pixel(src[x_stride] + offset_val);
                x_stride += stride;
            }
            width --;
        }

    }
    if (sao_eo_class != SAO_EO_HORIZ) {
        if (borders[1]){
            int offset_val = sao_offset_val[0];
            for (x = init_x; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + offset_val);
            }
            init_y = 1;
        }
        if (borders[3]){
            int offset_val = sao_offset_val[0];
            int y_stride   = stride * (_height-1);
            for (x = init_x; x < width; x++) {
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + offset_val);
            }
            height--;
        }
    }
    {
        int y_stride     = init_y * stride;
        int pos_0_0      = pos[sao_eo_class][0][0];
        int pos_0_1      = pos[sao_eo_class][0][1];
        int pos_1_0      = pos[sao_eo_class][1][0];
        int pos_1_1      = pos[sao_eo_class][1][1];

        int y_stride_0_1 = (init_y + pos_0_1) * stride;
        int y_stride_1_1 = (init_y + pos_1_1) * stride;
        for (y = init_y; y < height; y++) {
            for (x = init_x; x < width; x++) {
                int diff0         = CMP(src[x + y_stride], src[x + pos_0_0 + y_stride_0_1]);
                int diff1         = CMP(src[x + y_stride], src[x + pos_1_0 + y_stride_1_1]);
                int offset_val    = edge_idx[2 + diff0 + diff1];
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + sao_offset_val[offset_val]);
            }
            y_stride     += stride;
            y_stride_0_1 += stride;
            y_stride_1_1 += stride;
        }
    }
#undef CMP
}

static void FUNC(sao_edge_filter_wpp_1)(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx)
{
    int x, y;
    uint8_t *dst = _dst;   // put here pixel
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int chroma = c_idx!=0;
    //struct SAOParams *sao;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_eo_class = sao->eo_class[c_idx];

    const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };
    const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };

    int init_x = 0, init_y = 0, width = _width, height = _height;

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))

            init_y = -(4>>chroma)-2;
            if(!borders[2] )
                width -= ((8>>chroma)+2);
            height = (4>>chroma)+2;

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    init_y = init_x = 0;
    if (sao_eo_class != SAO_EO_VERT) {
        if (borders[0]) {
            int offset_val = sao_offset_val[0];
            int y_stride   = 0;
            for (y = 0; y < height; y++) {
                dst[y_stride] = av_clip_pixel(src[y_stride] + offset_val);
                y_stride += stride;
            }
            init_x = 1;
        }
        if (borders[2]) {
            int offset_val = sao_offset_val[0];
            int x_stride   = _width-1;
            for (x = 0; x < height; x++) {
                dst[x_stride] = av_clip_pixel(src[x_stride] + offset_val);
                x_stride += stride;
            }
            width --;
        }

    }

    {
        int y_stride     = init_y * stride;
        int pos_0_0      = pos[sao_eo_class][0][0];
        int pos_0_1      = pos[sao_eo_class][0][1];
        int pos_1_0      = pos[sao_eo_class][1][0];
        int pos_1_1      = pos[sao_eo_class][1][1];

        int y_stride_0_1 = (init_y + pos_0_1) * stride;
        int y_stride_1_1 = (init_y + pos_1_1) * stride;
        for (y = init_y; y < height; y++) {
            for (x = init_x; x < width; x++) {
                int diff0         = CMP(src[x + y_stride], src[x + pos_0_0 + y_stride_0_1]);
                int diff1         = CMP(src[x + y_stride], src[x + pos_1_0 + y_stride_1_1]);
                int offset_val    = edge_idx[2 + diff0 + diff1];
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + sao_offset_val[offset_val]);
            }
            y_stride     += stride;
            y_stride_0_1 += stride;
            y_stride_1_1 += stride;
        }
    }
#undef CMP
}

static void FUNC(sao_edge_filter_wpp_2)(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx)
{
    int x, y;
    uint8_t *dst = _dst;   // put here pixel
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int chroma = c_idx!=0;
    //struct SAOParams *sao;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_eo_class = sao->eo_class[c_idx];

    const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };
    const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };

    int init_x = 0, init_y = 0, width = _width, height = _height;

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))


            init_x = -(8>>chroma)-2;
            width = (8>>chroma)+2;
            if(!borders[3])
                height -= ((4>>chroma)+2);

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    init_y = init_x = 0;

    if (sao_eo_class != SAO_EO_HORIZ) {
        if (borders[1]){
            int offset_val = sao_offset_val[0];
            for (x = init_x; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + offset_val);
            }
            init_y = 1;
        }
        if (borders[3]){
            int offset_val = sao_offset_val[0];
            int y_stride   = stride * (_height-1);
            for (x = init_x; x < width; x++) {
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + offset_val);
            }
            height--;
        }
    }
    {
        int y_stride     = init_y * stride;
        int pos_0_0      = pos[sao_eo_class][0][0];
        int pos_0_1      = pos[sao_eo_class][0][1];
        int pos_1_0      = pos[sao_eo_class][1][0];
        int pos_1_1      = pos[sao_eo_class][1][1];

        int y_stride_0_1 = (init_y + pos_0_1) * stride;
        int y_stride_1_1 = (init_y + pos_1_1) * stride;
        for (y = init_y; y < height; y++) {
            for (x = init_x; x < width; x++) {
                int diff0         = CMP(src[x + y_stride], src[x + pos_0_0 + y_stride_0_1]);
                int diff1         = CMP(src[x + y_stride], src[x + pos_1_0 + y_stride_1_1]);
                int offset_val    = edge_idx[2 + diff0 + diff1];
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + sao_offset_val[offset_val]);
            }
            y_stride     += stride;
            y_stride_0_1 += stride;
            y_stride_1_1 += stride;
        }
    }
#undef CMP
}

static void FUNC(sao_edge_filter_wpp_3)(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, struct SAOParams *sao,int *borders, int _width, int _height, int c_idx)
{
    int x, y;
    uint8_t *dst = _dst;   // put here pixel
    uint8_t *src = _src;
    ptrdiff_t stride = _stride;
    int chroma = c_idx!=0;
    //struct SAOParams *sao;
    int *sao_offset_val = sao->offset_val[c_idx];
    int sao_eo_class = sao->eo_class[c_idx];

    const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };
    const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };

    int init_x = 0, init_y = 0, width = _width, height = _height;

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))


            init_y = -(4>>chroma)-2;
            init_x = -(8>>chroma)-2;
            width = (8>>chroma)+2;
            height = (4>>chroma)+2;

    dst = dst + (init_y*_stride + init_x);
    src = src + (init_y*_stride + init_x);
    init_y = init_x = 0;


    {
        int y_stride     = init_y * stride;
        int pos_0_0      = pos[sao_eo_class][0][0];
        int pos_0_1      = pos[sao_eo_class][0][1];
        int pos_1_0      = pos[sao_eo_class][1][0];
        int pos_1_1      = pos[sao_eo_class][1][1];

        int y_stride_0_1 = (init_y + pos_0_1) * stride;
        int y_stride_1_1 = (init_y + pos_1_1) * stride;
        for (y = init_y; y < height; y++) {
            for (x = init_x; x < width; x++) {
                int diff0         = CMP(src[x + y_stride], src[x + pos_0_0 + y_stride_0_1]);
                int diff1         = CMP(src[x + y_stride], src[x + pos_1_0 + y_stride_1_1]);
                int offset_val    = edge_idx[2 + diff0 + diff1];
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + sao_offset_val[offset_val]);
            }
            y_stride     += stride;
            y_stride_0_1 += stride;
            y_stride_1_1 += stride;
        }
    }
#undef CMP
}



static void FUNC(sao_band_filter)(uint8_t * _dst, uint8_t *_src, ptrdiff_t _stride, int *sao_offset_val,
                                  int sao_left_class, int width, int height)
{
    pixel *dst = (pixel*)_dst;
    pixel *src = (pixel*)_src;
    ptrdiff_t stride = _stride/sizeof(pixel);
    int band_table[32] = { 0 };
    int k, y, x;
    int shift = BIT_DEPTH - 5;

    for (k = 0; k < 4; k++)
        band_table[(k + sao_left_class) & 31] = k + 1;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
#ifndef GCC_OPTIMIZATION_ENABLE
            x++;
            dst[x] = av_clip_pixel(src[x] + sao_offset_val[band_table[src[x] >> shift]]);
#endif
        }
        dst += stride;
        src += stride;
    }
}

static void FUNC(sao_edge_filter)(uint8_t *_dst, uint8_t *_src, ptrdiff_t _stride, int *sao_offset_val,
                                  int sao_eo_class, int at_top_border, int at_bottom_border,
                                  int at_left_border, int at_right_border,
                                  int width, int height)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    pixel *src = (pixel*)_src;
    ptrdiff_t stride = _stride/sizeof(pixel);

    const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };
    const uint8_t edge_idx[] = { 1, 2, 0, 3, 4 };

    int init_x = 0, init_y = 0;

#ifndef OPTIMIZATION_ENABLE
    int border_edge_idx = 0;
#define DST(x, y) dst[(x) + stride * (y)]
#define SRC(x, y) src[(x) + stride * (y)]

#define FILTER(x, y, edge_idx)                                      \
    DST(x, y) = av_clip_pixel(SRC(x, y) + sao_offset_val[edge_idx])

#define DIFF(x, y, k) CMP(SRC(x, y), SRC((x) + pos[sao_eo_class][(k)][0],       \
                                         (y) + pos[sao_eo_class][(k)][1]))
#endif
#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))

    if (sao_eo_class != SAO_EO_VERT) {
        if (at_left_border) {
#ifdef OPTIMIZATION_ENABLE
            int offset_val = sao_offset_val[0];
            int y_stride   = 0;
            for (y = 0; y < height; y++) {
                dst[y_stride] = av_clip_pixel(src[y_stride] + offset_val);
                y_stride += stride;
            }
#else
            for (y = 0; y < height; y++)
                FILTER(0, y, border_edge_idx);
#endif
            init_x = 1;
        }
        if (at_right_border) {
#ifdef OPTIMIZATION_ENABLE
            int offset_val = sao_offset_val[0];
            int x_stride   = width - 1;
            for (x = 0; x < height; x++) {
                dst[x_stride] = av_clip_pixel(src[x_stride] + offset_val);
                x_stride += stride;
            }
#else
            for (x = 0; x < height; x++)
                FILTER(width - 1, x, border_edge_idx);
#endif
            width--;
        }
    }
    if (sao_eo_class != SAO_EO_HORIZ) {
        if (at_top_border) {
#ifdef OPTIMIZATION_ENABLE
            int offset_val = sao_offset_val[0];
            for (x = init_x; x < width; x++) {
                dst[x] = av_clip_pixel(src[x] + offset_val);
            }
#else
            for (x = init_x; x < width; x++)
                FILTER(x, 0, border_edge_idx);
#endif
            init_y = 1;
        }
        if (at_bottom_border) {
#ifdef OPTIMIZATION_ENABLE
            int offset_val = sao_offset_val[0];
            int y_stride   = stride * (height - 1);
            for (x = init_x; x < width; x++) {
                dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + offset_val);
            }
#else
            for (x = init_x; x < width; x++)
                FILTER(x, height - 1, border_edge_idx);
#endif
            height--;
        }
    }
#ifdef OPTIMIZATION_ENABLE
    {
    int y_stride     = init_y * stride;
    int pos_0_0      = pos[sao_eo_class][0][0];
    int pos_0_1      = pos[sao_eo_class][0][1];
    int pos_1_0      = pos[sao_eo_class][1][0];
    int pos_1_1      = pos[sao_eo_class][1][1];
    int y_stride_0_1 = (init_y + pos_0_1) * stride;
    int y_stride_1_1 = (init_y + pos_1_1) * stride;
    for (y = init_y; y < height; y++) {
        for (x = init_x; x < width; x++) {
            int diff0         = CMP(src[x + y_stride], src[x + pos_0_0 + y_stride_0_1]);
            int diff1         = CMP(src[x + y_stride], src[x + pos_1_0 + y_stride_1_1]);
            int offset_val    = edge_idx[2 + diff0 + diff1];
            dst[x + y_stride] = av_clip_pixel(src[x + y_stride] + sao_offset_val[offset_val]);
        }
        y_stride     += stride;
        y_stride_0_1 += stride;
        y_stride_1_1 += stride;
    }
    }
#else
    for (y = init_y; y < height; y++) {
        for (x = init_x; x < width; x++)
            FILTER(x, y, edge_idx[2 + DIFF(x, y, 0) + DIFF(x, y, 1)]);
    }
#endif

#ifndef OPTIMIZATION_ENABLE
#undef DST
#undef SRC
#undef FILTER
#undef DIFF
#endif
#undef CMP
}

#undef SET
#undef SCALE
#undef ADD_AND_SCALE
#undef TR_4
#undef TR_4_1
#undef TR_4_2
#undef TR_8
#undef TR_8_1
#undef TR_8_2
#undef TR_16
#undef TR_16_1
#undef TR_16_2
#undef TR_32
#undef TR_32_1
#undef TR_32_2


static void FUNC(put_hevc_qpel_pixels)(int16_t *dst, ptrdiff_t dststride,
                                       uint8_t *_src, ptrdiff_t _srcstride,
                                       int width, int height, int16_t* mcbuffer)
{
    int x, y;
    pixel *src = (pixel*)_src;
    ptrdiff_t srcstride = _srcstride/sizeof(pixel);
    
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = src[x] << (14 - BIT_DEPTH);
#ifndef GCC_OPTIMIZATION_ENABLE
            x++;
            dst[x] = src[x] << (14 - BIT_DEPTH);
            x++;
            dst[x] = src[x] << (14 - BIT_DEPTH);
            x++;
            dst[x] = src[x] << (14 - BIT_DEPTH);
#endif
        }
        src += srcstride;
        dst += dststride;
    }
}


#define QPEL_FILTER_1(src, stride)                                              \
    (-src[x-3*stride] + 4*src[x-2*stride] - 10*src[x-stride] + 58*src[x] +      \
     17*src[x+stride] - 5*src[x+2*stride] + 1*src[x+3*stride])
#define QPEL_FILTER_2(src, stride)                                              \
    (-src[x-3*stride] + 4*src[x-2*stride] - 11*src[x-stride] + 40*src[x] +      \
     40*src[x+stride] - 11*src[x+2*stride] + 4*src[x+3*stride] - src[x+4*stride])
#define QPEL_FILTER_3(src, stride)                                              \
    (src[x-2*stride] - 5*src[x-stride] + 17*src[x] + 58*src[x+stride]           \
     - 10*src[x+2*stride] + 4*src[x+3*stride] - src[x+4*stride])


#if  defined GCC_OPTIMIZATION_ENABLE
#define PUT_HEVC_QPEL_H(H)                                                      \
static void FUNC(put_hevc_qpel_h ## H)(int16_t *dst, ptrdiff_t dststride,       \
uint8_t *_src, ptrdiff_t _srcstride,  \
int width, int height, int16_t* mcbuffer)                \
{                                                                               \
int x, y;                                                                   \
pixel *src = (pixel*)_src;                                                  \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                             \
\
for (y = 0; y < height; y++) {                                              \
for (x = 0; x < width; x++)                                             \
dst[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);             \
src += srcstride;                                                       \
dst += dststride;                                                       \
}                                                                           \
}
#define PUT_HEVC_QPEL_V(V)                                                      \
static void FUNC(put_hevc_qpel_v ## V)(int16_t *dst, ptrdiff_t dststride,       \
uint8_t *_src, ptrdiff_t _srcstride,  \
int width, int height, int16_t mcbuffer)                \
{                                                                               \
int x, y;                                                                   \
pixel *src = (pixel*)_src;                                                  \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                             \
\
for (y = 0; y < height; y++)  {                                             \
for (x = 0; x < width; x++)                                             \
dst[x] = QPEL_FILTER_ ## V (src, srcstride) >> (BIT_DEPTH - 8);     \
src += srcstride;                                                       \
dst += dststride;                                                       \
}                                                                           \
}
#define PUT_HEVC_QPEL_HV(H, V)                                                            \
static void FUNC(put_hevc_qpel_h ## H ## v ## V )(int16_t *dst, ptrdiff_t dststride,      \
uint8_t *_src, ptrdiff_t _srcstride,    \
int width, int height, int16_t* mcbuffer)                  \
{                                                                                         \
int x, y;                                                                             \
pixel *src = (pixel*)_src;                                                            \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                                       \
\
/*int16_t tmp_array[(MAX_PB_SIZE+7)*MAX_PB_SIZE];  */                                     \
int16_t *tmp = mcbuffer;                                                             \
\
src -= qpel_extra_before[V] * srcstride;                                              \
\
for (y = 0; y < height + qpel_extra[V]; y++) {                                        \
for (x = 0; x < width; x++)                                                       \
tmp[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);                       \
src += srcstride;                                                                 \
tmp += MAX_PB_SIZE;                                                               \
}                                                                                     \
\
tmp = tmp_array + qpel_extra_before[V] * MAX_PB_SIZE;                                 \
\
for (y = 0; y < height; y++) {                                                        \
for (x = 0; x < width; x++)                                                       \
dst[x] = QPEL_FILTER_ ## V (tmp, MAX_PB_SIZE) >> 6;                           \
tmp += MAX_PB_SIZE;                                                               \
dst += dststride;                                                                 \
}                                                                                     \
}

#else
#define PUT_HEVC_QPEL_H(H)                                                      \
static void FUNC(put_hevc_qpel_h ## H)(int16_t *dst, ptrdiff_t dststride,       \
uint8_t *_src, ptrdiff_t _srcstride,  \
int width, int height, int16_t* mcbuffer)                \
{                                                                               \
int x, y;                                                                   \
pixel *src = (pixel*)_src;                                                  \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                             \
\
for (y = 0; y < height; y++) {                                              \
for (x = 0; x < width; x++) {                                           \
dst[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);             \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);             \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);             \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);             \
}                                                                       \
src += srcstride;                                                       \
dst += dststride;                                                       \
}                                                                           \
}
#define PUT_HEVC_QPEL_V(V)                                                      \
static void FUNC(put_hevc_qpel_v ## V)(int16_t *dst, ptrdiff_t dststride,       \
uint8_t *_src, ptrdiff_t _srcstride,  \
int width, int height, int16_t* mcbuffer)                \
{                                                                               \
int x, y;                                                                   \
pixel *src = (pixel*)_src;                                                  \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                             \
\
for (y = 0; y < height; y++)  {                                             \
for (x = 0; x < width; x++) {                                           \
dst[x] = QPEL_FILTER_ ## V (src, srcstride) >> (BIT_DEPTH - 8);     \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## V (src, srcstride) >> (BIT_DEPTH - 8);     \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## V (src, srcstride) >> (BIT_DEPTH - 8);     \
x++;                                                                \
dst[x] = QPEL_FILTER_ ## V (src, srcstride) >> (BIT_DEPTH - 8);     \
}                                                                       \
src += srcstride;                                                       \
dst += dststride;                                                       \
}                                                                           \
}
#define PUT_HEVC_QPEL_HV(H, V)                                                            \
static void FUNC(put_hevc_qpel_h ## H ## v ## V )(int16_t *dst, ptrdiff_t dststride,      \
uint8_t *_src, ptrdiff_t _srcstride,    \
int width, int height, int16_t* mcbuffer)                  \
{                                                                                         \
int x, y;                                                                             \
pixel *src = (pixel*)_src;                                                            \
ptrdiff_t srcstride = _srcstride/sizeof(pixel);                                       \
\
/*int16_t tmp_array[(MAX_PB_SIZE+7)*MAX_PB_SIZE];                                       */\
int16_t *tmp = mcbuffer;                                                             \
\
src -= qpel_extra_before[V] * srcstride;                                              \
\
for (y = 0; y < height + qpel_extra[V]; y++) {                                        \
for (x = 0; x < width; x++) {                                                     \
tmp[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);                       \
x++;                                                                          \
tmp[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);                       \
x++;                                                                          \
tmp[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);                       \
x++;                                                                          \
tmp[x] = QPEL_FILTER_ ## H (src, 1) >> (BIT_DEPTH - 8);                       \
}                                                                                 \
src += srcstride;                                                                 \
tmp += MAX_PB_SIZE;                                                               \
}                                                                                     \
\
tmp = mcbuffer + qpel_extra_before[V] * MAX_PB_SIZE;                                 \
\
for (y = 0; y < height; y++) {                                                        \
for (x = 0; x < width; x++) {                                                     \
dst[x] = QPEL_FILTER_ ## V (tmp, MAX_PB_SIZE) >> 6;                           \
x++;                                                                          \
dst[x] = QPEL_FILTER_ ## V (tmp, MAX_PB_SIZE) >> 6;                           \
x++;                                                                          \
dst[x] = QPEL_FILTER_ ## V (tmp, MAX_PB_SIZE) >> 6;                           \
x++;                                                                          \
dst[x] = QPEL_FILTER_ ## V (tmp, MAX_PB_SIZE) >> 6;                           \
}                                                                                 \
tmp += MAX_PB_SIZE;                                                               \
dst += dststride;                                                                 \
}                                                                                     \
}
#endif

PUT_HEVC_QPEL_H(1)
PUT_HEVC_QPEL_H(2)
PUT_HEVC_QPEL_H(3)
PUT_HEVC_QPEL_V(1)
PUT_HEVC_QPEL_V(2)
PUT_HEVC_QPEL_V(3)
PUT_HEVC_QPEL_HV(1, 1)
PUT_HEVC_QPEL_HV(1, 2)
PUT_HEVC_QPEL_HV(1, 3)
PUT_HEVC_QPEL_HV(2, 1)
PUT_HEVC_QPEL_HV(2, 2)
PUT_HEVC_QPEL_HV(2, 3)
PUT_HEVC_QPEL_HV(3, 1)
PUT_HEVC_QPEL_HV(3, 2)
PUT_HEVC_QPEL_HV(3, 3)

static void FUNC(put_hevc_epel_pixels)(int16_t *dst, ptrdiff_t dststride,
                                       uint8_t *_src, ptrdiff_t _srcstride,
                                       int width, int height, int mx, int my, int16_t* mcbuffer)
{
    int x, y;
    pixel *src = (pixel*)_src;
    ptrdiff_t srcstride = _srcstride/sizeof(pixel);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = src[x] << (14 - BIT_DEPTH);
#ifndef GCC_OPTIMIZATION_ENABLE
            x++;
            dst[x] = src[x] << (14 - BIT_DEPTH);
#endif
        }
        src += srcstride;
        dst += dststride;
    }
}


#ifdef OPTIMIZATION_ENABLE
#define EPEL_FILTER(src, stride) \
    (filter_0*src[x-stride] + filter_1*src[x] + filter_2*src[x+stride] + filter_3*src[x+2*stride])
#else
#define EPEL_FILTER(src, stride, F) \
    (F[0]*src[x-stride] + F[1]*src[x] + F[2]*src[x+stride] + F[3]*src[x+2*stride])
#endif

static void FUNC(put_hevc_epel_h)(int16_t *dst, ptrdiff_t dststride,
                                  uint8_t *_src, ptrdiff_t _srcstride,
                                  int width, int height, int mx, int my, int16_t* mcbuffer)
{
    int x, y;
    pixel *src = (pixel*)_src;
    ptrdiff_t srcstride = _srcstride/sizeof(pixel);
    const int8_t *filter = epel_filters[mx-1];
#ifdef OPTIMIZATION_ENABLE
    int8_t filter_0 = filter[0];
    int8_t filter_1 = filter[1];
    int8_t filter_2 = filter[2];
    int8_t filter_3 = filter[3];
#endif
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
#ifdef GCC_OPTIMIZATION_ENABLE
            dst[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
#else
            dst[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
            x++;
            dst[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
#endif
        }
        src += srcstride;
        dst += dststride;
    }
}

static void FUNC(put_hevc_epel_v)(int16_t *dst, ptrdiff_t dststride,
                                  uint8_t *_src, ptrdiff_t _srcstride,
                                  int width, int height, int mx, int my, int16_t* mcbuffer)
{
    int x, y;
    pixel *src = (pixel*)_src;
    ptrdiff_t srcstride = _srcstride/sizeof(pixel);
    const int8_t *filter = epel_filters[my-1];
#ifdef OPTIMIZATION_ENABLE
    int8_t filter_0 = filter[0];
    int8_t filter_1 = filter[1];
    int8_t filter_2 = filter[2];
    int8_t filter_3 = filter[3];
#endif
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
#ifdef GCC_OPTIMIZATION_ENABLE
            dst[x] = EPEL_FILTER(src, srcstride) >> (BIT_DEPTH - 8);
#else
            dst[x] = EPEL_FILTER(src, srcstride) >> (BIT_DEPTH - 8);
            x++;
            dst[x] = EPEL_FILTER(src, srcstride) >> (BIT_DEPTH - 8);
#endif
        }
        src += srcstride;
        dst += dststride;
    }
}

static void FUNC(put_hevc_epel_hv)(int16_t *dst, ptrdiff_t dststride,
                                   uint8_t *_src, ptrdiff_t _srcstride,
                                   int width, int height, int mx, int my, int16_t* mcbuffer)
{
    int x, y;
    pixel *src = (pixel*)_src;
    ptrdiff_t srcstride = _srcstride/sizeof(pixel);
    const int8_t *filter_h = epel_filters[mx-1];
    const int8_t *filter_v = epel_filters[my-1];
#ifdef OPTIMIZATION_ENABLE
    int8_t filter_0 = filter_h[0];
    int8_t filter_1 = filter_h[1];
    int8_t filter_2 = filter_h[2];
    int8_t filter_3 = filter_h[3];
#endif
    //#ifndef OPTIMIZATION_ENABLE
   // int16_t tmp_array[(MAX_PB_SIZE+3)*MAX_PB_SIZE];
    //#endif
    int16_t *tmp = mcbuffer;
    
    src -= epel_extra_before * srcstride;

    for (y = 0; y < height + epel_extra; y++) {
        for (x = 0; x < width; x++) {
#ifdef GCC_OPTIMIZATION_ENABLE
            tmp[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
#else
            tmp[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
            x++;
            tmp[x] = EPEL_FILTER(src, 1) >> (BIT_DEPTH - 8);
#endif
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }
    
    tmp = mcbuffer + epel_extra_before * MAX_PB_SIZE;
#ifdef OPTIMIZATION_ENABLE
    filter_0 = filter_v[0];
    filter_1 = filter_v[1];
    filter_2 = filter_v[2];
    filter_3 = filter_v[3];
#endif
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
#ifdef GCC_OPTIMIZATION_ENABLE
            dst[x] = EPEL_FILTER(tmp, MAX_PB_SIZE) >> 6;
#else
            dst[x] = EPEL_FILTER(tmp, MAX_PB_SIZE) >> 6;
            x++;
            dst[x] = EPEL_FILTER(tmp, MAX_PB_SIZE) >> 6;
#endif
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void FUNC(put_unweighted_pred)(uint8_t *_dst, ptrdiff_t _dststride,
                                      int16_t *src, ptrdiff_t srcstride,
                                      int width, int height)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t dststride = _dststride/sizeof(pixel);

    int shift = 14 - BIT_DEPTH;
#if BIT_DEPTH < 14
    int offset = 1 << (shift - 1);
#else
    int offset = 0;
#endif
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = av_clip_pixel((src[x] + offset) >> shift);
#ifndef GCC_OPTIMIZATION_ENABLE
            x++;
            dst[x] = av_clip_pixel((src[x] + offset) >> shift);
            x++;
            dst[x] = av_clip_pixel((src[x] + offset) >> shift);
            x++;
            dst[x] = av_clip_pixel((src[x] + offset) >> shift);
#endif
        }
        dst += dststride;
        src += srcstride;
    }
}


static void FUNC(put_weighted_pred_avg)(uint8_t *_dst, ptrdiff_t _dststride,
                                        int16_t *src1, int16_t *src2, ptrdiff_t srcstride,
                                        int width, int height)
{
    int x, y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t dststride = _dststride/sizeof(pixel);

    int shift = 14 + 1 - BIT_DEPTH;
#if BIT_DEPTH < 14
    int offset = 1 << (shift - 1);
#else
    int offset = 0;
#endif
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = av_clip_pixel((src1[x] + src2[x] + offset) >> shift);
#ifndef GCC_OPTIMIZATION_ENABLE
            x++;
            dst[x] = av_clip_pixel((src1[x] + src2[x] + offset) >> shift);
            x++;
            dst[x] = av_clip_pixel((src1[x] + src2[x] + offset) >> shift);
            x++;
            dst[x] = av_clip_pixel((src1[x] + src2[x] + offset) >> shift);
#endif
        }
        dst  += dststride;
        src1 += srcstride;
        src2 += srcstride;
    }
}


static void FUNC(weighted_pred)(uint8_t denom, int16_t wlxFlag, int16_t olxFlag,
                                     uint8_t *_dst, ptrdiff_t _dststride,
                                     int16_t *src, ptrdiff_t srcstride,
                                     int width, int height)
{
    int shift;
    int log2Wd;
    int wx;
    int ox;
    int x , y;
    int offset;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t dststride = _dststride/sizeof(pixel);

    shift = 14 - BIT_DEPTH;
    log2Wd = denom + shift;
    offset = 1 << (log2Wd - 1);
    wx = wlxFlag;
    ox = olxFlag * ( 1 << ( BIT_DEPTH - 8 ) );

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            if (log2Wd >= 1) {
                dst[x] = av_clip_pixel(((src[x] * wx + offset) >> log2Wd) + ox);
            } else {
                dst[x] = av_clip_pixel(src[x] * wx + ox);
            }
        }
        dst  += dststride;
        src  += srcstride;
    }
}

static void FUNC(weighted_pred_avg)(uint8_t denom, int16_t wl0Flag, int16_t wl1Flag,
                                         int16_t ol0Flag, int16_t ol1Flag, uint8_t *_dst, ptrdiff_t _dststride,
                                         int16_t *src1, int16_t *src2, ptrdiff_t srcstride,
                                         int width, int height)
{
    int shift;
    int log2Wd;
    int w0;
    int w1;
    int o0;
    int o1;
    int x , y;
    pixel *dst = (pixel*)_dst;
    ptrdiff_t dststride = _dststride/sizeof(pixel);

    shift = 14 - BIT_DEPTH;
    log2Wd = denom + shift;
    w0 = wl0Flag;
    w1 = wl1Flag;
    o0 = (ol0Flag) * ( 1 << ( BIT_DEPTH - 8 ) );
    o1 = (ol1Flag) * ( 1 << ( BIT_DEPTH - 8 ) );

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = av_clip_pixel((src1[x] * w0 + src2[x] * w1 + ((o0 + o1 + 1) << log2Wd)) >> (log2Wd + 1));
        }
        dst  += dststride;
        src1 += srcstride;
        src2 += srcstride;
    }
}


// line zero
#define P3 pix[-4*xstride]
#define P2 pix[-3*xstride]
#define P1 pix[-2*xstride]
#define P0 pix[-xstride]
#define Q0 pix[0]
#define Q1 pix[xstride]
#define Q2 pix[2*xstride]
#define Q3 pix[3*xstride]

// line three. used only for deblocking decision
#define TP3 pix[-4*xstride+3*ystride]
#define TP2 pix[-3*xstride+3*ystride]
#define TP1 pix[-2*xstride+3*ystride]
#define TP0 pix[-xstride+3*ystride]
#define TQ0 pix[3*ystride]
#define TQ1 pix[xstride+3*ystride]
#define TQ2 pix[2*xstride+3*ystride]
#define TQ3 pix[3*xstride+3*ystride]

static void FUNC(hevc_loop_filter_luma_h)(uint8_t *_pix, ptrdiff_t _xstride, ptrdiff_t _ystride,
                                        int no_p, int no_q, int _beta, int _tc)
{
    int d;
    pixel *pix = (pixel*)_pix;
    ptrdiff_t xstride = _xstride/sizeof(pixel);
    ptrdiff_t ystride = _ystride/sizeof(pixel);
    const int dp0 = abs(P2 - 2 * P1 +  P0);
    const int dq0 = abs(Q2 - 2 * Q1 +  Q0);
    const int dp3 = abs(TP2 - 2 * TP1 + TP0);
    const int dq3 = abs(TQ2 - 2 * TQ1 + TQ0);
    const int d0 = dp0 + dq0;
    const int d3 = dp3 + dq3;

    const int beta = _beta << (BIT_DEPTH - 8);
    const int tc = _tc << (BIT_DEPTH - 8);

    if (d0 + d3 < beta) {
        const int beta_3 = beta >> 3;
        const int beta_2 = beta >> 2;
        const int tc25 = ((tc * 5 + 1) >> 1);

        if(abs( P3 -  P0) + abs( Q3 -  Q0) < beta_3 && abs( P0 -  Q0) < tc25 &&
           abs(TP3 - TP0) + abs(TQ3 - TQ0) < beta_3 && abs(TP0 - TQ0) < tc25 &&
                                 (d0 << 1) < beta_2 &&      (d3 << 1) < beta_2) {
            // strong filtering
            const int tc2 = tc << 1;
            for(d = 0; d < 4; d++) {
                const int p3 = P3;
                const int p2 = P2;
                const int p1 = P1;
                const int p0 = P0;
                const int q0 = Q0;
                const int q1 = Q1;
                const int q2 = Q2;
                const int q3 = Q3;
                if(!no_p) {
                    P0 = av_clip_c(( p2 + 2*p1 + 2*p0 + 2*q0 + q1 + 4 ) >> 3, p0-tc2, p0+tc2);
                    P1 = av_clip_c(( p2 + p1 + p0 + q0 + 2 ) >> 2, p1-tc2, p1+tc2);
                    P2 = av_clip_c(( 2*p3 + 3*p2 + p1 + p0 + q0 + 4 ) >> 3, p2-tc2, p2+tc2);
                }
                if(!no_q) {
                    Q0 = av_clip_c(( p1 + 2*p0 + 2*q0 + 2*q1 + q2 + 4 ) >> 3, q0-tc2, q0+tc2);
                    Q1 = av_clip_c(( p0 + q0 + q1 + q2 + 2 ) >> 2, q1-tc2, q1+tc2);
                    Q2 = av_clip_c(( 2*q3 + 3*q2 + q1 + q0 + p0 + 4 ) >> 3, q2-tc2, q2+tc2);
                }
                pix += ystride;
            }
        } else { // normal filtering
            int nd_p = 1;
            int nd_q = 1;
            const int tc_2 = tc >> 1;
            if (dp0 + dp3 < ((beta+(beta>>1))>>3))
                nd_p = 2;
            if (dq0 + dq3 < ((beta+(beta>>1))>>3))
                nd_q = 2;

            for(d = 0; d < 4; d++) {
                const int p2 = P2;
                const int p1 = P1;
                const int p0 = P0;
                const int q0 = Q0;
                const int q1 = Q1;
                const int q2 = Q2;
                int delta0 = (9*(q0 - p0) - 3*(q1 - p1) + 8) >> 4;
                if (abs(delta0) < 10 * tc) {
                    delta0 = av_clip_c(delta0, -tc, tc);
                    if(!no_p)
                        P0 = av_clip_pixel(p0 + delta0);
                    if(!no_q)
                        Q0 = av_clip_pixel(q0 - delta0);
                    if(!no_p && nd_p > 1) {
                        const int deltap1 = av_clip_c((((p2 + p0 + 1) >> 1) - p1 + delta0) >> 1, -tc_2, tc_2);
                        P1 = av_clip_pixel(p1 + deltap1);
                    }
                    if(!no_q && nd_q > 1) {
                        const int deltaq1 = av_clip_c((((q2 + q0 + 1) >> 1) - q1 - delta0) >> 1, -tc_2, tc_2);
                        Q1 = av_clip_pixel(q1 + deltaq1);
                    }
                }
                pix += ystride;
            }
        }
    }
}

static void FUNC(hevc_loop_filter_luma_v)(uint8_t *_pix, ptrdiff_t _xstride, ptrdiff_t _ystride,
                                        int no_p, int no_q, int _beta, int _tc)
{
    int d;
    pixel *pix = (pixel*)_pix;
    ptrdiff_t xstride = _xstride/sizeof(pixel);
    ptrdiff_t ystride = _ystride/sizeof(pixel);
    const int dp0 = abs(P2 - 2 * P1 +  P0);
    const int dq0 = abs(Q2 - 2 * Q1 +  Q0);
    const int dp3 = abs(TP2 - 2 * TP1 + TP0);
    const int dq3 = abs(TQ2 - 2 * TQ1 + TQ0);
    const int d0 = dp0 + dq0;
    const int d3 = dp3 + dq3;

    const int beta = _beta << (BIT_DEPTH - 8);
    const int tc = _tc << (BIT_DEPTH - 8);

    if (d0 + d3 < beta) {
        const int beta_3 = beta >> 3;
        const int beta_2 = beta >> 2;
        const int tc25 = ((tc * 5 + 1) >> 1);

        if(abs( P3 -  P0) + abs( Q3 -  Q0) < beta_3 && abs( P0 -  Q0) < tc25 &&
           abs(TP3 - TP0) + abs(TQ3 - TQ0) < beta_3 && abs(TP0 - TQ0) < tc25 &&
                                 (d0 << 1) < beta_2 &&      (d3 << 1) < beta_2) {
            // strong filtering
            const int tc2 = tc << 1;
            for(d = 0; d < 4; d++) {
                const int p3 = P3;
                const int p2 = P2;
                const int p1 = P1;
                const int p0 = P0;
                const int q0 = Q0;
                const int q1 = Q1;
                const int q2 = Q2;
                const int q3 = Q3;
                if(!no_p) {
                    P0 = av_clip_c(( p2 + 2*p1 + 2*p0 + 2*q0 + q1 + 4 ) >> 3, p0-tc2, p0+tc2);
                    P1 = av_clip_c(( p2 + p1 + p0 + q0 + 2 ) >> 2, p1-tc2, p1+tc2);
                    P2 = av_clip_c(( 2*p3 + 3*p2 + p1 + p0 + q0 + 4 ) >> 3, p2-tc2, p2+tc2);
                }
                if(!no_q) {
                    Q0 = av_clip_c(( p1 + 2*p0 + 2*q0 + 2*q1 + q2 + 4 ) >> 3, q0-tc2, q0+tc2);
                    Q1 = av_clip_c(( p0 + q0 + q1 + q2 + 2 ) >> 2, q1-tc2, q1+tc2);
                    Q2 = av_clip_c(( 2*q3 + 3*q2 + q1 + q0 + p0 + 4 ) >> 3, q2-tc2, q2+tc2);
                }
                pix += ystride;
            }
        } else { // normal filtering
            int nd_p = 1;
            int nd_q = 1;
            const int tc_2 = tc >> 1;
            if (dp0 + dp3 < ((beta+(beta>>1))>>3))
                nd_p = 2;
            if (dq0 + dq3 < ((beta+(beta>>1))>>3))
                nd_q = 2;

            for(d = 0; d < 4; d++) {
                const int p2 = P2;
                const int p1 = P1;
                const int p0 = P0;
                const int q0 = Q0;
                const int q1 = Q1;
                const int q2 = Q2;
                int delta0 = (9*(q0 - p0) - 3*(q1 - p1) + 8) >> 4;
                if (abs(delta0) < 10 * tc) {
                    delta0 = av_clip_c(delta0, -tc, tc);
                    if(!no_p)
                        P0 = av_clip_pixel(p0 + delta0);
                    if(!no_q)
                        Q0 = av_clip_pixel(q0 - delta0);
                    if(!no_p && nd_p > 1) {
                        const int deltap1 = av_clip_c((((p2 + p0 + 1) >> 1) - p1 + delta0) >> 1, -tc_2, tc_2);
                        P1 = av_clip_pixel(p1 + deltap1);
                    }
                    if(!no_q && nd_q > 1) {
                        const int deltaq1 = av_clip_c((((q2 + q0 + 1) >> 1) - q1 - delta0) >> 1, -tc_2, tc_2);
                        Q1 = av_clip_pixel(q1 + deltaq1);
                    }
                }
                pix += ystride;
            }
        }
    }
}

static void FUNC(hevc_loop_filter_chroma)(uint8_t *_pix, ptrdiff_t _xstride, ptrdiff_t _ystride, int no_p, int no_q, int _tc)
{
    int d;
    pixel *pix = (pixel*)_pix;
    ptrdiff_t xstride = _xstride/sizeof(pixel);
    ptrdiff_t ystride = _ystride/sizeof(pixel);
    const int tc = _tc << (BIT_DEPTH - 8);

    for(d = 0; d < 4; d++) {
        int delta0;
        const int p1 = P1;
        const int p0 = P0;
        const int q0 = Q0;
        const int q1 = Q1;
        delta0 = av_clip_c((((q0 - p0) << 2) + p1 - q1 + 4) >> 3, -tc, tc);
        if(!no_p)
            P0 = av_clip_pixel(p0 + delta0);
        if(!no_q)
            Q0 = av_clip_pixel(q0 - delta0);
        pix += ystride;
    }
}

#undef P3
#undef P2
#undef P1
#undef P0
#undef Q0
#undef Q1
#undef Q2
#undef Q3

#undef TP3
#undef TP2
#undef TP1
#undef TP0
#undef TQ0
#undef TQ1
#undef TQ2
#undef TQ3







#ifdef SVC_EXTENSION
#define LumHor_FILTER(pel, coeff) \
(pel[0]*coeff[0] + pel[1]*coeff[1] + pel[2]*coeff[2] + pel[3]*coeff[3] + pel[4]*coeff[4] + pel[5]*coeff[5] + pel[6]*coeff[6] + pel[7]*coeff[7])

#define CroHor_FILTER(pel, coeff) \
(pel[0]*coeff[0] + pel[1]*coeff[1] + pel[2]*coeff[2] + pel[3]*coeff[3])

#define LumVer_FILTER(pel, coeff, stride) \
(pel[0]*coeff[0] + pel[stride]*coeff[1] + pel[2*stride]*coeff[2] + pel[3*stride]*coeff[3] + pel[4*stride]*coeff[4] + pel[5*stride]*coeff[5] + pel[6*stride]*coeff[6] + pel[7*stride]*coeff[7])

#define CroVer_FILTER(pel, coeff, stride) \
(pel[0]*coeff[0] + pel[stride]*coeff[1] + pel[2*stride]*coeff[2] + pel[3*stride]*coeff[3])
// Define the function for up-sampling


static void FUNC(upsample_base_picture)(AVFrame *FrameEL, AVFrame *FrameBL, AVFrame *Buffer, const int32_t enabled_up_sample_filter_luma[16][8], const int32_t enabled_up_sample_filter_chroma[12][4], struct HEVCWindow *Enhscal)
{
    int i,j;
    
#if SCALED_REF_LAYER_OFFSETS
    int widthBL =  FrameBL->width;
    int heightBL = FrameBL->height;
    int strideBL = FrameBL->linesize[0];
    int widthEL =  FrameEL->width - Enhscal->left_offset - Enhscal->right_offset;
    int heightEL = FrameEL->height - Enhscal->top_offset - Enhscal->bottom_offset;
    int strideEL = FrameEL->linesize[0];
#else
    const HEVCWindow &confBL; //= pcBasePic->getConformanceWindow();
    const HEVCWindow &confEL; //= pcUsPic->getConformanceWindow();
    
    int widthBL   = FrameBL->width - confBL.left_offset - confBL.right_offset;
    int heightBL  = FrameBL->height - confBL.top_offset - confBL.bottom_offset;
    int strideBL  = FrameBL->linesize[0];
    
    int widthEL   = FrameEL->width - confEL.getWindowLeftOffset() - confEL.getWindowRightOffset();
    int heightEL  = FrameEL->height - confEL.getWindowTopOffset() - confEL.getWindowBottomOffset();
    int strideEL  = FrameEL->linesize[0];
#endif
    pixel *srcBufY = (pixel*)FrameBL->data[0];;
    pixel *dstBufY = (pixel*)FrameEL->data[0];
    pixel *tempBufY = (pixel*)Buffer->data[0];
    pixel *srcY, *dstY;
    
    pixel *srcBufU = (pixel*)FrameBL->data[1];
    pixel *dstBufU = (pixel*)FrameEL->data[1];
    pixel *tempBufU = (pixel*)Buffer->data[1];
    pixel *srcU, *dstU;
    
    pixel *srcBufV = (pixel*)FrameBL->data[2];
    pixel *dstBufV = (pixel*)FrameEL->data[2];
    pixel *tempBufV = (pixel*)Buffer->data[2];
    pixel *srcV, *dstV;
    
    if( widthEL == widthBL && heightEL == heightBL ) //Copy frame
    {
        srcY = srcBufY - Enhscal->left_offset - Enhscal->right_offset * strideEL;
        dstY = dstBufY;
        for( i = 0; i < heightEL; i++ )
        {
            memcpy( dstY, srcY, sizeof(pixel) * widthBL );
            srcY += strideBL;
            srcY += strideEL;
        }
        
        widthEL  >>= 1;
        heightEL >>= 1;
        
        widthBL  >>= 1;
        heightBL >>= 1;
        
        strideBL  = FrameBL->linesize[1];
        strideEL  = FrameEL->linesize[1];
        
        srcU = srcBufU - ( Enhscal->left_offset >> 1 ) - ( Enhscal->top_offset >> 1 ) * strideEL;
        srcV = srcBufV - ( Enhscal->left_offset >> 1 ) - ( Enhscal->top_offset >> 1 ) * strideEL;
        
        dstU = dstBufU;
        dstV = dstBufV;
        
        for( i = 0; i < heightEL; i++ )
        {
            memcpy( dstU, srcU, sizeof(pixel) * widthBL );
            memcpy( dstV, srcV, sizeof(pixel) * widthBL );
            srcU += strideBL;
            srcV += strideBL;
            dstU += strideEL;
            dstV += strideEL;
        }
    }
    else
    {
        printf("######## Upsample the base layer %d %d %d %d %d %d \n", widthBL, heightBL, strideBL, widthEL, heightEL, strideEL);
#if PHASE_DERIVATION_IN_INTEGER
        int refPos16 = 0;
        int phase    = 0;
        int refPos   = 0;
        int32_t* coeff = enabled_up_sample_filter_chroma[phase];
#else
        int chromaPhaseDenominator;
        if (widthEL == 2*widthBL) // 2x scalability
            chromaPhaseDenominator = 8;
        else
            chromaPhaseDenominator = 12;
#endif
        
        
        
        
#if PHASE_DERIVATION_IN_INTEGER
        int   shiftX = 16;
        int   shiftY = 16;
        
        int   phaseX = 0;
        int   phaseY = 0;
        
        int   addX       = ( ( ( widthBL * phaseX ) << ( shiftX - 2 ) ) + ( widthEL >> 1 ) ) / widthEL + ( 1 << ( shiftX - 5 ) );
        int   addY       = ( ( ( heightBL * phaseY ) << ( shiftY - 2 ) ) + ( heightEL >> 1 ) ) / heightEL+ ( 1 << ( shiftY - 5 ) );
        
        int   deltaX     = 4 * phaseX;
        int   deltaY     = 4 * phaseY;
        
        int shiftXM4 = shiftX - 4;
        int shiftYM4 = shiftY - 4;
        
        int   scaleX     = ( ( widthBL << shiftX ) + ( widthEL >> 1 ) ) / widthEL;
        int   scaleY     = ( ( heightBL << shiftY ) + ( heightEL >> 1 ) ) / heightEL;
#else
        const double sFactor = 1.0 * widthBL / widthEL;
        const double sFactor12 = sFactor * 12;
#endif
#if ILP_DECODED_PICTURE
        widthEL   = FrameEL->width;  //pcUsPic->getWidth ();
        heightEL  = FrameEL->height; //pcUsPic->getHeight();
        
        widthBL   = FrameBL->width;
        heightBL  = FrameBL->height <= heightEL ? FrameBL->height:heightEL;  // min( FrameBL->height, heightEL);
#endif
#if SCALED_REF_LAYER_OFFSETS
        int leftStartL = Enhscal->left_offset;
        int rightEndL  = FrameEL->width - Enhscal->right_offset;
        int topStartL  = Enhscal->top_offset;
        int bottomEndL = FrameEL->height - Enhscal->bottom_offset;
#endif
        
        
        
        for( i = 0; i < widthEL; i++ )
        {
#if SCALED_REF_LAYER_OFFSETS
            int x = av_clip_c(i, leftStartL, rightEndL);
#endif
#if PHASE_DERIVATION_IN_INTEGER
#if SCALED_REF_LAYER_OFFSETS
            refPos16 = (((x - leftStartL)*scaleX + addX) >> shiftXM4) - deltaX;
#else
            refPos16 = ((i*scaleX + addX) >> shiftXM4) - deltaX;
#endif
            phase    = refPos16 & 15;
            refPos   = refPos16 >> 4;
            coeff = enabled_up_sample_filter_luma[phase];
#else
            int refPos12 = (Int) ( i * sFactor12 );
            int refPos = (Int)( i * sFactor );
            int phase = (refPos12 + 12) % 12;
            int32_t* coeff = enabled_up_sample_filter_luma[phase];
#endif
            
            srcY = srcBufY + refPos -((NTAPS_LUMA>>1) - 1);
            dstY = tempBufY + i;
            
            for( j = 0; j < heightBL ; j++ )
            {
                *dstY = LumHor_FILTER(srcY, coeff);
                srcY += strideBL;
                dstY += strideEL;
            }
        }
        
        
        //========== vertical upsampling ===========
        /*TempPic->setBorderExtension(false);
         pcTempPic->setHeight(heightBL);
         pcTempPic->extendPicBorder   (); // extend the border.
         pcTempPic->setHeight(heightEL);*/
        
        Buffer->height = heightEL;
        Buffer->width = heightBL;
        
        
        const int nShift = US_FILTER_PREC*2;
        int iOffset = 1 << (nShift - 1);
        
#if SCALED_REF_LAYER_OFFSETS
        for( j = 0; j < Buffer->height; j++ )
#else
            for( j = 0; j < heightEL; j++ )
#endif
            {
#if SCALED_REF_LAYER_OFFSETS
                int y = av_clip_c(j, topStartL, bottomEndL-1);
#endif
#if PHASE_DERIVATION_IN_INTEGER
#if SCALED_REF_LAYER_OFFSETS
                refPos16 = ((( y - topStartL )*scaleY + addY) >> shiftYM4) - deltaY;
#else
                refPos16 = ((j*scaleY + addY) >> shiftYM4) - deltaY;
#endif
                phase    = refPos16 & 15;
                refPos   = refPos16 >> 4;
                coeff = enabled_up_sample_filter_luma[phase];
#else
                Int refPos12 = (Int) (j * sFactor12 );
                Int refPos = (Int)( j * sFactor );
                Int phase = (refPos12 + 12) % 12;
                Int* coeff = enabled_up_sample_filter_luma[phase];
#endif
                
                srcY = tempBufY + (refPos -((NTAPS_LUMA>>1) - 1))*strideEL;
                dstY = dstBufY + j * strideEL;
#if SCALED_REF_LAYER_OFFSETS
                for( i = 0; i < Buffer->width; i++ )
#else
                    for( i = 0; i < widthEL; i++ )
#endif
                    {
                        *dstY = av_clip_pixel( (LumVer_FILTER(srcY, coeff, strideEL) + iOffset) >> (nShift));
#if SCALED_REF_LAYER_OFFSETS
                        // Only increase the x position of reference upsample picture when within the window
                        // "-2" to ensure that pointer doesn't go beyond the boundary rightEndL-1
                        if( (i >= leftStartL) && (i <= rightEndL-2) )
                        {
                            srcY++;
                        }
#else
                        srcY++;
#endif
                        dstY++;
                    }
            }
        
#if ILP_DECODED_PICTURE
#if SCALED_REF_LAYER_OFFSETS
        widthBL   = FrameBL->width; //pcBasePic->getWidth ();
        heightBL  = FrameBL->height; // pcBasePic->getHeight();
        
        widthEL   = FrameEL->width - Enhscal->right_offset - Enhscal->left_offset;
        heightEL  = FrameEL->height - Enhscal->top_offset - Enhscal->bottom_offset; //pcUsPic->getHeight() - scalEL.getWindowTopOffset()  - scalEL.getWindowBottomOffset();
#else
        widthBL   = FrameBL->width; //- confBL.getWindowLeftOffset() - confBL.getWindowRightOffset();
        heightBL  = FrameEL->height; //- confBL.getWindowTopOffset() - confBL.getWindowBottomOffset();
        
        widthEL   = FrameEL->width; //- confEL.getWindowLeftOffset() - confEL.getWindowRightOffset();
        heightEL  = FrameEL->height; //- confEL.getWindowTopOffset() - confEL.getWindowBottomOffset();
#endif
#endif
        
        //========== UV component upsampling ===========
        
        widthEL  >>= 1;
        heightEL >>= 1;
        
        widthBL  >>= 1;
        heightBL >>= 1;
        
        strideBL  = FrameBL->linesize[1];
        strideEL  = FrameEL->linesize[1];
#if SCALED_REF_LAYER_OFFSETS
        int leftStartC = Enhscal->left_offset>>1;
        int rightEndC  = (FrameEL->width>>1) - (Enhscal->right_offset>>1);
        int topStartC  = Enhscal->top_offset>>1;
        int bottomEndC = (FrameEL->height>>1) - (Enhscal->bottom_offset>>1);
#endif
        
#if PHASE_DERIVATION_IN_INTEGER
        shiftX = 16;
        shiftY = 16;
        
        phaseX = 0;
        phaseY = 1;
        
        addX       = ( ( ( widthBL * phaseX ) << ( shiftX - 2 ) ) + ( widthEL >> 1 ) ) / widthEL + ( 1 << ( shiftX - 5 ) );
        addY       = ( ( ( heightBL * phaseY ) << ( shiftY - 2 ) ) + ( heightEL >> 1 ) ) / heightEL+ ( 1 << ( shiftY - 5 ) );
        
        deltaX     = 4 * phaseX;
        deltaY     = 4 * phaseY;
        
        shiftXM4 = shiftX - 4;
        shiftYM4 = shiftY - 4;
        
        scaleX     = ( ( widthBL << shiftX ) + ( widthEL >> 1 ) ) / widthEL;
        scaleY     = ( ( heightBL << shiftY ) + ( heightEL >> 1 ) ) / heightEL;
#endif
        
#if ILP_DECODED_PICTURE
        widthEL   = FrameEL->width >> 1;
        heightEL  = FrameEL->height >> 1;
        
        widthBL   = FrameBL->width >> 1;
        //        heightBL  = min( FrameBL->height, heightEL );
        heightBL  = FrameBL->height <= heightEL ? FrameBL->height:heightEL;
#endif
        
        //========== horizontal upsampling ===========
        for( i = 0; i < widthEL; i++ )
        {
#if SCALED_REF_LAYER_OFFSETS
            int x = av_clip_c(i, leftStartC, rightEndC - 1);
#endif
#if PHASE_DERIVATION_IN_INTEGER
#if SCALED_REF_LAYER_OFFSETS
            refPos16 = (((x - leftStartC)*scaleX + addX) >> shiftXM4) - deltaX;
#else
            refPos16 = ((i*scaleX + addX) >> shiftXM4) - deltaX;
#endif
            phase    = refPos16 & 15;
            refPos   = refPos16 >> 4;
            coeff = enabled_up_sample_filter_chroma[phase];
#else
            int refPosM = (Int) ( i * chromaPhaseDenominator * sFactor );
            int refPos = (Int)( i * sFactor );
            int phase = (refPosM + chromaPhaseDenominator) % chromaPhaseDenominator;
            int* coeff = enabled_up_sample_filter_chroma[phase];
#endif
            
            srcU = srcBufU + refPos -((NTAPS_CHROMA>>1) - 1);
            srcV = srcBufV + refPos -((NTAPS_CHROMA>>1) - 1);
            dstU = tempBufU + i;
            dstV = tempBufV + i;
            
            for( j = 0; j < heightBL ; j++ )
            {
                *dstU = CroHor_FILTER(srcU, coeff);
                *dstV = CroHor_FILTER(srcV, coeff);
                
                srcU += strideBL;
                srcV += strideBL;
                dstU += strideEL;
                dstV += strideEL;
            }
        }
        
        //========== vertical upsampling ===========
        /* pcTempPic->setBorderExtension(false);
         pcTempPic->setHeight(heightBL << 1);
         pcTempPic->extendPicBorder   (); // extend the border.
         pcTempPic->setHeight(heightEL << 1);*/
        
#if SCALED_REF_LAYER_OFFSETS
        for( j = 0; j < Buffer->height >> 1; j++ )
#else
            for( j = 0; j < heightEL; j++ )
#endif
            {
#if SCALED_REF_LAYER_OFFSETS
                int y = av_clip_c(j, topStartC, bottomEndC - 1);
#endif
#if PHASE_DERIVATION_IN_INTEGER
#if SCALED_REF_LAYER_OFFSETS
                refPos16 = (((y - topStartC)*scaleY + addY) >> shiftYM4) - deltaY;
#else
                refPos16 = ((j*scaleY + addY) >> shiftYM4) - deltaY;
#endif
                phase    = refPos16 & 15;
                refPos   = refPos16 >> 4;
                coeff = enabled_up_sample_filter_chroma[phase];
#else
                Int refPosM = (Int) (j * chromaPhaseDenominator * sFactor) - 1;
                Int refPos;
                if ( refPosM < 0 )
                {
                    refPos = (Int)( j * sFactor ) - 1;
                }
                else
                {
                    refPos = refPosM / chromaPhaseDenominator;
                }
                int phase = (refPosM + chromaPhaseDenominator) % chromaPhaseDenominator;
                int* coeff = enabled_up_sample_filter_chroma[phase];
#endif
                
                srcU = tempBufU  + (refPos -((NTAPS_CHROMA>>1) - 1))*strideEL;
                srcV = tempBufV  + (refPos -((NTAPS_CHROMA>>1) - 1))*strideEL;
                
                dstU = dstBufU + j*strideEL;
                dstV = dstBufV + j*strideEL;
                
#if SCALED_REF_LAYER_OFFSETS
                for( i = 0; i < Buffer->width >> 1; i++ )
#else
                    for( i = 0; i < widthEL; i++ )
#endif
                    {
                        *dstU = av_clip_pixel( (CroVer_FILTER(srcU, coeff, strideEL) + iOffset) >> (nShift));
                        *dstV = av_clip_pixel( (CroVer_FILTER(srcV, coeff, strideEL) + iOffset) >> (nShift));
                        
                        
#if SCALED_REF_LAYER_OFFSETS
                        // Only increase the x position of reference upsample picture when within the window
                        // "-2" to ensure that pointer doesn't go beyond the boundary rightEndC-1
                        if( (i >= leftStartC) && (i <= rightEndC-2) )
                        {
                            srcU++;
                            srcV++;
                        }
#else
                        srcU++;
                        srcV++;
#endif
                        dstU++;
                        dstV++;
                    }
            }
    }
}


#undef LumHor_FILTER
#undef LumCro_FILTER
#undef LumVer_FILTER
#undef CroVer_FILTER
#endif

