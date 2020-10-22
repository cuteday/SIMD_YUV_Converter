#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <mmintrin.h>
#include <emmintrin.h>
#include <x86intrin.h>


using namespace std;

#define DEFAULT_WIDTH   1920
#define DEFAULT_HEIGHT  1080

// Channel synatx
#define YUV_Y   0
#define YUV_U   1
#define YUV_V   2

#define ARGB_R   0
#define ARGB_G   1
#define ARGB_B   2
#define ARGB_A   3

typedef unsigned char byte_t;

typedef union {
    short i16[4];
    __m64 pkg;
} m64;

typedef union {
    short i16[8];
    float f32[4];
    __m128 pkg;
} m128;

typedef union {
    short i16[16];
    float f32[8];
    __m256 pkg;
} m256;

#define NUM_MODES   1
typedef enum{
    SIMD_NULL,
    SIMD_MMX,
    SIMD_SSE2,
    SIMD_AVX,
} SIMD_MODE;

class Graph{
    // Photo meta-data
public:
    byte_t *graph;
    int width, height;
    int fsize, npixel;
    inline byte_t &operator[](int index) { return graph[index]; }
    virtual inline byte_t& operator() (int r, int c, int channel) = 0;
};

class ARGB: public Graph{
public:
    ARGB(int w = DEFAULT_WIDTH, int h = DEFAULT_HEIGHT){
        width = w, height = h;
        fsize = w * h * 4;
        graph = (byte_t*)calloc(fsize, sizeof(byte_t));
    }
    inline byte_t& operator() (int r, int c, int channel){
        return graph[r * width * 4 + c * 4 + channel];
    }
};

class YUV: public Graph{
public:
    int nuv;        // number of bytes of U / V takes
    byte_t *graph;
    YUV(int w = DEFAULT_WIDTH, int h = DEFAULT_HEIGHT){
        width = w, height = h;
        nuv = (w / 2) * (h / 2), npixel = w * h;
        fsize = w * h + 2 * nuv;
        graph = (byte_t *)calloc(fsize, sizeof(byte_t));
    }
    void load(const char *path){
        // fprintf(stdout, "Loading file from %s...\n", path);
        std::ifstream file(path, std::ios::binary | std::ios::in);
        file.read((char *)graph, fsize);
        file.close();
    }
    void save(const char *path){
        // for (int i = 0;i<fsize;i++)
        //     printf("%d ", graph[i]);
        fprintf(stdout, "Saving YUV file to %s...\n", path);
        std::ofstream file(path, std::ios::binary | std::ios::out);
        file.write((char *)graph, fsize);
        file.close();
    }
    inline byte_t& operator() (int r, int c, int channel){
        if(channel == YUV_Y)
            return graph[r * width + c];
        if(channel == YUV_U)
            return graph[npixel + (r / 2) * (width / 2) + c / 2];
        return graph[npixel + nuv + (r / 2) * (width / 2) + c / 2];
    }
};

/* ___________________________ Convertion Functions ___________________________ */

extern void YUV_ARGB(YUV& yuv, ARGB& argb, byte_t alpha, SIMD_MODE mode);
extern void ALPHA_MIX(ARGB& argb, SIMD_MODE mode);
extern void RGB_YUV(ARGB& argb, YUV& yuv, SIMD_MODE mode);

/* ___________________________ Parameters List ___________________________ */

#define PR_RY   1.164383
#define PR_RU   0.
#define PR_RV   1.596027
#define PR_GY   1.164383
#define PR_GU  -0.391762
#define PR_GV  -0.812968
#define PR_BY   1.164383
#define PR_BU   2.017232
#define PR_BV   0.
#define CO_R    -222.921584
#define CO_G    135.575312
#define CO_B    -276.835824

#define PR_YR   0.256788
#define PR_YG   0.504129
#define PR_YB   0.097906
#define PR_UR  -0.148223
#define PR_UG  -0.290993
#define PR_UB   0.439216
#define PR_VR   0.439216
#define PR_VG  -0.367788
#define PR_VB  -0.071427
#define CO_Y    16.
#define CO_U    128.
#define CO_V    128.

extern float param_yuv2rgb[3][4];
extern float param_rgb2yuv[3][4];
extern const char *simd_names[];
