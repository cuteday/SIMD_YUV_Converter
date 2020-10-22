#include "function.hpp"

float param_yuv2rgb[3][4] = {
    {PR_RY, PR_RU, PR_RV, CO_R },
    {PR_GY, PR_GU, PR_GV, CO_G },
    {PR_BY, PR_BU, PR_BV, CO_B },
};

float param_rgb2yuv[3][4] = {
    {PR_YR, PR_YG, PR_YB, CO_Y},
    {PR_UR, PR_UG, PR_UB, CO_U},
    {PR_VR, PR_VG, PR_VB, CO_V},
};

const char *simd_names[] = {
    "None SIMD",
    "MMX",
    "SSE2",
    "AVX",
};