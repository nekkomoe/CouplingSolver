// 将物理场插值或下采样

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <functional>
using namespace std;

void withfile(string filename, function<void(FILE*)> f) {
    FILE *fp = fopen(filename.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: cannot open file %s\n", filename.c_str());
        exit(1);
    }
    f(fp);
    fclose(fp);
}

void downsample(int n, int m, double *src,
            int tar_n, int tar_m, double *dst) {
    // src[n][m] -> dst[tar_n][tar_m]
    for (int i = 0; i < tar_n; i++) {
        for (int j = 0; j < tar_m; j++) {
            int x = i * n / tar_n;
            int y = j * m / tar_m;
            dst[i*tar_m+j] = src[x*m+y];
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <data_path>\n", argv[0]);
        return 0;
    }
    string data_path = argv[1];
    string file_pf = data_path + "/pf_final.bin";
    string file_lc = data_path + "/lc_final.bin";
    int Nr, Nz;
    int new_Nr = 32, new_Nz = 32;
    withfile(data_path+"/info.txt", [&](FILE *fp) {
        fscanf(fp, "%d%d", &Nr, &Nz);
    });
    int scale[5][2] = {
        {Nr, Nz},   // Ts
        {Nr, Nz},   // Tc
        {Nr, Nz},   // Phi
        {Nr, Nz},   // P
        {Nr, Nz+1}, // U
    };
    int scale_lamC[5][2] = {
        {Nr, Nz},   // lambda*C
    };
    withfile(file_pf, [&](FILE *fp) {
        double pf[Nr*Nz*4+Nr*(Nz+1)];
        fread(pf, sizeof(double), Nr*Nz*4+Nr*(Nz+1), fp);
        double new_pf[new_Nr*new_Nz*4+new_Nr*(new_Nz+1)];
        for (int i = 0; i < 5; i++) {
            int n = scale[i][0], m = scale[i][1];
            int new_n = new_Nr, new_m = new_Nz;
            if (i == 4) new_m++;
            downsample(n, m, pf+i*Nr*Nz, new_n, new_m, new_pf+i*new_Nr*new_Nz);
        }
        FILE *fp_new = fopen((data_path+"/pf_final_downsampled.bin").c_str(), "wb");
        fwrite(new_pf, sizeof(double), new_Nr*new_Nz*4+new_Nr*(new_Nz+1), fp_new);
        fclose(fp_new);
    });
    withfile(file_lc, [&](FILE *fp) {
        double lc[Nr*Nz];
        fread(lc, sizeof(double), Nr*Nz, fp);
        double new_lc[new_Nr*new_Nz];
        downsample(Nr, Nz, lc, new_Nr, new_Nz, new_lc);
        FILE *fp_new = fopen((data_path+"/lc_final_downsampled.bin").c_str(), "wb");
        fwrite(new_lc, sizeof(double), new_Nr*new_Nz, fp_new);
        fclose(fp_new);
    });
}