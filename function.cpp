#include "function.hpp"
using namespace std;
// #define LOOP_UNROLLING
// #define LIMIT(x) (max(min(x, 255.), 0.))
inline byte_t LIMIT(short x) { return max((short)0, min((short)255, x)); }

void YUV_ARGB(YUV& yuv, ARGB& argb, byte_t alpha, SIMD_MODE mode){
    assert(0 <= alpha && alpha <= 255);
    int w = yuv.width, h = yuv.height;
    if(mode==SIMD_NULL){
        for (int i = 0; i < h;i++)
            for (int j = 0; j < w;j++){
                // argb(i, j, 0) = LIMIT(param_yuv2rgb[0][0] * yuv(i, j, 0) + param_yuv2rgb[0][1] * yuv(i, j, 1) + param_yuv2rgb[0][2] * yuv(i, j, 2) + param_yuv2rgb[0][3]);
                // argb(i, j, 1) = LIMIT(param_yuv2rgb[1][0] * yuv(i, j, 0) + param_yuv2rgb[1][1] * yuv(i, j, 1) + param_yuv2rgb[1][2] * yuv(i, j, 2) + param_yuv2rgb[1][3]);
                // argb(i, j, 2) = LIMIT(param_yuv2rgb[2][0] * yuv(i, j, 0) + param_yuv2rgb[2][1] * yuv(i, j, 1) + param_yuv2rgb[2][2] * yuv(i, j, 2) + param_yuv2rgb[2][3]);
                for (int ch = 0; ch < 4;ch++)
                    argb(i, j, ch) = LIMIT(param_yuv2rgb[ch][0] * yuv(i, j, 0) + param_yuv2rgb[ch][1] * yuv(i, j, 1) + param_yuv2rgb[ch][2] * yuv(i, j, 2) + param_yuv2rgb[ch][3]);
                argb(i, j, 3) = alpha;
            }
    }
    else if(mode==SIMD_MMX){
        short param_mmx[3][4];     // 近似放大后的转换系数，用矩阵表示
        for (int i = 0; i < 3;i++)
            for (int j = 0; j < 3;j++)
                param_mmx[i][j] = (1 << 8) * param_yuv2rgb[i][j];
            
#ifdef LOOP_UNROLLING
        __m64 mult_y[3], mult_u[3], mult_v[3], constant[3];
        

        for (int i = 0; i < 3;i++){     // R G B 的 Y U V
            mult_y[i] = _mm_setr_pi16(param_mmx[i][0], param_mmx[i][0], param_mmx[i][0], param_mmx[i][0]);
            mult_u[i] = _mm_setr_pi16(param_mmx[i][1], param_mmx[i][1], param_mmx[i][1], param_mmx[i][1]);
            mult_v[i] = _mm_setr_pi16(param_mmx[i][2], param_mmx[i][2], param_mmx[i][2], param_mmx[i][2]);
            constant[i] = _mm_setr_pi16(param_yuv2rgb[i][3],param_yuv2rgb[i][3],param_yuv2rgb[i][3],param_yuv2rgb[i][3]);
        }

        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j += 4){
                ushort Y[] = {yuv(i, j, 0), yuv(i, j + 1, 0), yuv(i, j + 2, 0), yuv(i, j + 3, 0)};
                ushort U[] = {yuv(i, j, 1), yuv(i, j + 1, 1), yuv(i, j + 2, 1), yuv(i, j + 3, 1)};
                ushort V[] = {yuv(i, j, 2), yuv(i, j + 1, 2), yuv(i, j + 2, 2), yuv(i, j + 3, 2)};

                __m64 pkg_y = _mm_setr_pi16(Y[0], Y[1], Y[2], Y[3]);
                __m64 pkg_u = _mm_setr_pi16(U[0], U[1], U[2], U[3]);
                __m64 pkg_v = _mm_setr_pi16(V[0], V[1], V[2], V[3]);

                argb(i, j + 1, 3) = argb(i, j + 2, 3) = argb(i, j + 3, 3) = argb(i, j + 0, 3) = alpha;

                for (int ch = 0; ch < 3; ch++){     // R,G,B 3 Channels
                    __m64 res_y = _mm_mullo_pi16(pkg_y, mult_y[ch]);
                    __m64 res_u = _mm_mullo_pi16(pkg_u, mult_u[ch]);
                    __m64 res_v = _mm_mullo_pi16(pkg_v, mult_v[ch]);

                    __m64 res = _mm_add_pi16(res_y, res_u);
                    res = _mm_add_pi16(res, res_v);
                    res = _mm_srli_pi16(res, 8);
                    res = _mm_add_pi16(res, constant[ch]);

                    argb(i, j, ch) = (((m64 *)&res)->i16[0]), argb(i, j+1, ch) = (((m64 *)&res)->i16[1]),
                    argb(i, j+2, ch) = (((m64 *)&res)->i16[2]), argb(i, j+3, ch) = (((m64 *)&res)->i16[3]);
                }
            }
#else
        __m64 mult_yuv[3][3], constant[3];
        for (int i = 0; i < 3;i++){     // R G B loop
            constant[i] = _mm_setr_pi16(param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3]);
            for (int j = 0; j < 3; j++) // Y U V loop
                mult_yuv[i][j] = _mm_setr_pi16(param_mmx[i][j], param_mmx[i][j], param_mmx[i][j], param_mmx[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j+=4){
                for (int ch = 0; ch < 3;ch++){      // RGB loop
                    __m64 res = _mm_set_pi16(0, 0, 0, 0);
                    for (int k = 0; k < 3; k++){        // YUV loop
                        byte_t yuv_val[4] = { yuv(i, j, k), yuv(i, j + 1, k), yuv(i, j + 2, k), yuv(i, j + 3, k)};
                        __m64 yuv_pkg = _mm_setr_pi16(yuv_val[0], yuv_val[1], yuv_val[2], yuv_val[3]);
                        __m64 cur_res = _mm_mullo_pi16(mult_yuv[ch][k], yuv_pkg);
                        res = _mm_add_pi16(res, cur_res);
                    }
                    res = _mm_srai_pi16(res, 8);
                    res = _mm_add_pi16(res, constant[ch]);
                    argb(i, j, ch) = (((m64 *)&res)->i16[0]), argb(i, j+1, ch) = (((m64 *)&res)->i16[1]),
                    argb(i, j+2, ch) = (((m64 *)&res)->i16[2]), argb(i, j+3, ch) = (((m64 *)&res)->i16[3]);    
                }
                argb(i, j + 1, 3) = argb(i, j + 2, 3) = argb(i, j + 3, 3) = argb(i, j + 0, 3) = alpha;
            }

#endif
    }
    else if(mode==SIMD_SSE2){
        __m128 mult_yuv[3][3], constant[3];
        for (int i = 0; i < 3;i++){     // R G B loop
            constant[i] = _mm_set_ps(param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3]);
            for (int j = 0; j < 3; j++) // Y U V loop
                mult_yuv[i][j] = _mm_set_ps(param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j+=4){
                for (int ch = 0; ch < 3;ch++){      // RGB loop
                    __m128 res = _mm_set_ps(0., 0., 0., 0.);
                    for (int k = 0; k < 3; k++){        // YUV loop
                        byte_t yuv_val[4] = { yuv(i, j, k), yuv(i, j + 1, k), yuv(i, j + 2, k), yuv(i, j + 3, k)};
                        __m128 yuv_pkg = _mm_setr_ps((float)yuv_val[0], (float)yuv_val[1], (float)yuv_val[2], (float)yuv_val[3]);
                        __m128 cur_res = _mm_mul_ps(mult_yuv[ch][k], yuv_pkg);
                        res = _mm_add_ps(res, cur_res);
                    }
                    res = _mm_add_ps(res, constant[ch]);
                    argb(i, j, ch) = LIMIT(((m128*)&res)->f32[0]), argb(i, j+1, ch) = LIMIT(((m128 *)&res)->f32[1]),
                    argb(i, j+2, ch) = LIMIT(((m128 *)&res)->f32[2]), argb(i, j+3, ch) = LIMIT(((m128 *)&res)->f32[3]);    
                }
                argb(i, j + 1, 3) = argb(i, j + 2, 3) = argb(i, j + 3, 3) = argb(i, j + 0, 3) = alpha;
            }
    }
    else if(mode==SIMD_AVX){
        __m256 mult_yuv[3][3], constant[3];
        for (int i = 0; i < 3;i++){     // R G B loop
            constant[i] = _mm256_set_ps(param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3],
                                    param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3], param_yuv2rgb[i][3]);
            for (int j = 0; j < 3; j++) // Y U V loop
                mult_yuv[i][j] = _mm256_set_ps(param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j],
                                            param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j], param_yuv2rgb[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j+=8){
                for (int ch = 0; ch < 3;ch++){      // RGB loop
                    __m256 res = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
                    for (int k = 0; k < 3; k++){        // YUV loop
                        byte_t yuv_val[8] = { yuv(i, j, k), yuv(i, j + 1, k), yuv(i, j + 2, k), yuv(i, j + 3, k), 
                                            yuv(i, j + 4, k), yuv(i, j + 5, k), yuv(i, j + 6, k), yuv(i, j + 7, k)};
                        __m256 yuv_pkg = _mm256_setr_ps((float)yuv_val[0], (float)yuv_val[1], (float)yuv_val[2], (float)yuv_val[3],
                                                        (float)yuv_val[4], (float)yuv_val[5], (float)yuv_val[6], (float)yuv_val[7]);
                        __m256 cur_res = _mm256_mul_ps(mult_yuv[ch][k], yuv_pkg);
                        res = _mm256_add_ps(res, cur_res);
                    }
                    res = _mm256_add_ps(res, constant[ch]);
                    for (int k = 0; k < 8;k++)
                        argb(i, j + k, ch) = LIMIT(((m256 *)&res)->f32[k]);
                }
                for (int k = 0; k < 8;k++)
                        argb(i, j + k, 3) = alpha;
            } 
    }
    
}

void ALPHA_MIX(ARGB& argb, SIMD_MODE mode){
    int w = argb.width, h = argb.height;
    if(mode==SIMD_NULL){
        for (int i = 0; i < h;i++)
            for (int j = 0; j < w;j++)
                for (int ch = 0; ch < 3;ch++){
                    argb(i, j, ch) *= (argb(i, j, 3) / 255.);
                }
    }
    else if(mode==SIMD_MMX){
        __m64 pkg[3], res[3];
        for (int i = 0; i < argb.fsize;i+=16){
            __m64 alpha = _mm_setr_pi16(argb[i + 3], argb[i + 7], argb[i + 11], argb[i + 15]);
            for (int j = 0; j < 3; j++){
                pkg[j] = _mm_setr_pi16(argb[i + j + 0], argb[i + j + 4], argb[i + j + 8], argb[i + j + 12]);
                res[j] = _mm_srli_pi16(_mm_mullo_pi16(pkg[j], alpha), 8);
                argb[i + j + 0] = ((m64 *)&res[j])->i16[0] ,  argb[i + j + 4] = ((m64 *)&res[j])->i16[1],
                argb[i + j + 8] = ((m64 *)&res[j])->i16[2],  argb[i + j + 12] = ((m64 *)&res[j])->i16[3];
            }
        }
    }
    else if(mode==SIMD_SSE2){
        __m128i pkg, res;
        for (int i = 0; i < argb.fsize; i += 32){
            __m128i alpha = _mm_setr_epi16(argb[i + 3], argb[i + 7], argb[i + 11], argb[i + 15],
                                          argb[i + 19], argb[i + 23], argb[i + 27], argb[i + 31]);
            for (int j = 0; j < 3;j++){
                pkg = _mm_setr_epi16(argb[i + j + 0],argb[i + j + 4],argb[i + j + 8],argb[i + j + 12],
                                    argb[i + j + 16],argb[i + j + 20],argb[i + j + 24],argb[i + j + 28]);
                res = _mm_srli_epi16(_mm_mullo_epi16(pkg, alpha), 8);
                argb[i + j + 0] = ((m128 *)&res)->i16[0] ,  argb[i + j + 4] = ((m128 *)&res)->i16[1],
                argb[i + j + 8] = ((m128 *)&res)->i16[2],  argb[i + j + 12] = ((m128 *)&res)->i16[3],
                argb[i + j + 16] = ((m128 *)&res)->i16[4] ,  argb[i + j + 20] = ((m128 *)&res)->i16[5],
                argb[i + j + 24] = ((m128 *)&res)->i16[6],  argb[i + j + 28] = ((m128 *)&res)->i16[7];
            }
        }
    }
    else if(mode==SIMD_AVX){
        float factor = 1.0 / 256.;
        __m256 div = _mm256_set_ps(factor, factor, factor, factor, factor, factor, factor, factor);
        __m256 pkg, res;
        for (int i = 0; i < argb.fsize; i += 32){
            __m256 alpha = _mm256_setr_ps(argb[i + 3], argb[i + 7], argb[i + 11], argb[i + 15],
                                           argb[i + 19], argb[i + 23], argb[i + 27], argb[i + 31]);
            for (int j = 0; j < 3; j++){
                pkg = _mm256_setr_ps(argb[i + j + 0], argb[i + j + 4], argb[i + j + 8], argb[i + j + 12],
                                     argb[i + j + 16], argb[i + j + 20], argb[i + j + 24], argb[i + j + 28]);
                res = _mm256_mul_ps(_mm256_mul_ps(pkg, alpha), div);
                for (int k = 0; k < 8; k++)
                    argb[i + j + 4 * k] = ((m256 *)&res)->f32[k];
            }
        } 
    }
}

void RGB_YUV(ARGB& argb, YUV& yuv, SIMD_MODE mode){
    int w = argb.width, h = argb.height;
    if(mode==SIMD_NULL){
        for (int i = 0; i < h;i++)
            for (int j = 0; j < w;j++){
                // yuv(i,j,0) = LIMIT(param_rgb2yuv[0][0] * argb(i, j, 0) + param_rgb2yuv[0][1] * argb(i, j, 1) + param_rgb2yuv[0][2] * argb(i, j, 2) + param_rgb2yuv[0][3]);
                // yuv(i,j,1) = LIMIT(param_rgb2yuv[1][0] * argb(i, j, 0) + param_rgb2yuv[1][1] * argb(i, j, 1) + param_rgb2yuv[1][2] * argb(i, j, 2) + param_rgb2yuv[1][3]);
                // yuv(i,j,2) = LIMIT(param_rgb2yuv[2][0] * argb(i, j, 0) + param_rgb2yuv[2][1] * argb(i, j, 1) + param_rgb2yuv[2][2] * argb(i, j, 2) + param_rgb2yuv[2][3]);
                for (int ch = 0; ch < 3;ch++)
                    yuv(i,j,ch) = LIMIT(param_rgb2yuv[ch][0] * argb(i, j, 0) + param_rgb2yuv[ch][1] * argb(i, j, 1) + param_rgb2yuv[ch][2] * argb(i, j, 2) + param_rgb2yuv[ch][3]);
            }
    }
    else if(mode==SIMD_MMX){
        short param_mmx[3][4];     // 近似放大后的转换系数，用矩阵表示
        for (int i = 0; i < 3;i++)
            for (int j = 0; j < 3;j++)
                param_mmx[i][j] = (1 << 8) * param_rgb2yuv[i][j];

        __m64 mult_rgb[3][3], constant[3];
        for (int i = 0; i < 3;i++){ // YUV loop
            constant[i] = _mm_setr_pi16(param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3]);
            for (int j = 0; j < 3; j++) // RGB loop
                mult_rgb[i][j] = _mm_setr_pi16(param_mmx[i][j], param_mmx[i][j], param_mmx[i][j], param_mmx[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j += 4){
                for (int ch = 0; ch < 3; ch++){ // YUV loop
                    __m64 res = _mm_set_pi16(0, 0, 0, 0);
                    for (int k = 0; k < 3; k++){ // RGB loop
                        byte_t rgb_val[4] = {argb(i, j, k), argb(i, j + 1, k), argb(i, j + 2, k), argb(i, j + 3, k)};
                        __m64 rgb_pkg = _mm_setr_pi16(rgb_val[0], rgb_val[1], rgb_val[2], rgb_val[3]);
                        __m64 cur_res = _mm_mullo_pi16(mult_rgb[ch][k], rgb_pkg);
                        res = _mm_add_pi16(res, cur_res);
                    }
                    res = _mm_srli_pi16(res, 8);
                    res = _mm_adds_pi16(res, constant[ch]);

                    yuv(i, j, ch) = (((m64 *)&res)->i16[0]), yuv(i, j + 1, ch) = (((m64 *)&res)->i16[1]),
                    yuv(i, j + 2, ch) = (((m64 *)&res)->i16[2]), yuv(i, j + 3, ch) = (((m64 *)&res)->i16[3]);
                }
            }
    }
    else if(mode==SIMD_SSE2){
        __m128 mult_rgb[3][3], constant[3];
        for (int i = 0; i < 3;i++){ // YUV loop
            constant[i] = _mm_set_ps(param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3]);
            for (int j = 0; j < 3; j++) // RGB loop
                mult_rgb[i][j] = _mm_set_ps(param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j += 4){
                for (int ch = 0; ch < 3; ch++){ // YUV loop
                    __m128 res = _mm_set_ps(0., 0., 0., 0.);
                    for (int k = 0; k < 3; k++){ // RGB loop
                        byte_t rgb_val[4] = {argb(i, j, k), argb(i, j + 1, k), argb(i, j + 2, k), argb(i, j + 3, k)};
                        __m128 rgb_pkg = _mm_setr_ps((float)rgb_val[0], (float)rgb_val[1], (float)rgb_val[2], (float)rgb_val[3]);
                        __m128 cur_res = _mm_mul_ps(mult_rgb[ch][k], rgb_pkg);
                        res = _mm_add_ps(res, cur_res);
                    }
                    res = _mm_add_ps(res, constant[ch]);
                    yuv(i, j, ch) = LIMIT(((m128 *)&res)->f32[0]), yuv(i, j + 1, ch) = LIMIT(((m128 *)&res)->f32[1]),
                    yuv(i, j + 2, ch) = LIMIT(((m128 *)&res)->f32[2]), yuv(i, j + 3, ch) = LIMIT(((m128 *)&res)->f32[3]);
                }
            }
    }
    else if(mode==SIMD_AVX){
        __m256 mult_rgb[3][3], constant[3];
        for (int i = 0; i < 3;i++){     // R G B loop
            constant[i] = _mm256_set_ps(param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3],
                                    param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3], param_rgb2yuv[i][3]);
            for (int j = 0; j < 3; j++) // Y U V loop
                mult_rgb[i][j] = _mm256_set_ps(param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j],
                                            param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j], param_rgb2yuv[i][j]);
        }
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j+=8){
                for (int ch = 0; ch < 3;ch++){      // RGB loop
                    __m256 res = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
                    for (int k = 0; k < 3; k++){        // YUV loop
                        byte_t rgb_val[8] = { argb(i, j, k), argb(i, j + 1, k), argb(i, j + 2, k), argb(i, j + 3, k), 
                                            argb(i, j+ 4, k), argb(i, j + 5, k), argb(i, j + 6, k), argb(i, j + 7, k)};
                        __m256 rgb_pkg = _mm256_setr_ps((float)rgb_val[0], (float)rgb_val[1], (float)rgb_val[2], (float)rgb_val[3],
                                                        (float)rgb_val[4], (float)rgb_val[5], (float)rgb_val[6], (float)rgb_val[7]);
                        __m256 cur_res = _mm256_mul_ps(mult_rgb[ch][k], rgb_pkg);
                        res = _mm256_add_ps(res, cur_res);
                    }
                    res = _mm256_add_ps(res, constant[ch]);
                    for (int k = 0; k < 8; k++)
                        yuv(i, j + k, ch) = LIMIT(((m256 *)&res)->f32[k]);
                }
            } 
    }
}
