#include "function.hpp"
#include <ctime>
using namespace std;

const char *input_file = "./yuv/dem1.yuv";
const char *output_file = "./output.yuv";
int width = DEFAULT_WIDTH, height = DEFAULT_HEIGHT;
byte_t alpha_value = 255;
SIMD_MODE simd_mode = SIMD_NULL;

void run_task(const char* in, const char *out, byte_t alpha, SIMD_MODE mode){
    YUV input = YUV(width, height), result = YUV(width, height);
    ARGB argb = ARGB(width, height);
    input.load(in);
    YUV_ARGB(input, argb, alpha, mode);
    ALPHA_MIX(argb, mode);
    RGB_YUV(argb, result, mode);
    result.save(out);
}

void run_test(const char* test_name, SIMD_MODE smode);

int main(int argc, char* argv[]){
    for (int i = 0; i < argc;i++){
        if(!strcmp(argv[i], "-f"))
            input_file = argv[++i];
        if(!strcmp(argv[i], "-o"))
            output_file = argv[++i];
        if(!strcmp(argv[i], "-a"))
            alpha_value = atoi(argv[++i]);
        if(!strcmp(argv[i], "-s")){
            width = atoi(argv[++i]);
            height = atoi(argv[++i]);
        }
        if(!strcmp(argv[i], "-mmx"))
            simd_mode = SIMD_MMX;
        if(!strcmp(argv[i], "-sse2"))
            simd_mode = SIMD_SSE2;
        if(!strcmp(argv[i], "-avx"))
            simd_mode = SIMD_AVX;           
        if (!strcmp(argv[i], "-test")){        
            run_test("AVX", SIMD_AVX);
            // run_test("SSE2", SIMD_SSE2);
            // run_test("MMX", SIMD_MMX);
            // run_test("NONE", SIMD_NULL);
            exit(0);
        }
    }
    fprintf(stdout, "Starting convertion! Currently under %s mode\n", simd_names[simd_mode]);
    run_task(input_file, output_file, alpha_value, simd_mode);
    fprintf(stdout, "Convertion Succeed! Quitting...\n");
}

void run_test(const char* test_name, SIMD_MODE smode){
    int start = clock();
    fprintf(stdout, "Starting SIMD test! Test mode %s\n", test_name);
    for (int alpha = 1; alpha < 255; alpha += 3){
        char *file_name = new char[20];
        sprintf(file_name, "./output/%s_alpha%d.yuv", test_name, alpha);
        run_task(input_file, file_name, alpha, smode);
    }
    fprintf(stdout, "%s Test finished! File saved to /output.\nTotal Time(sec): %.02f", test_name, (clock() - start) / 1e6);
}