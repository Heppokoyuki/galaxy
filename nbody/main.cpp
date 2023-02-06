#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <omp.h>
#include <sstream>
#include <chrono>
#include <immintrin.h>

// #define DEBUG_TIME

const float eps = 0.01;
const float eps2 = eps*eps;
const float MAX_T = 10;

using namespace std;

union particle {
    float f[8];
    struct  {
        float x, y, z, m, vx, vy, vz, eps2;
    } __attribute__((aligned(32)));
};

string
getFileNameWithoutExt(string f)
{
    int ext_i = f.find_last_of(".");
    return f.substr(0, ext_i);
}

string
createFileName(string input_file)
{
    /* create time string */
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    char time_string[64] = "";
    sprintf(time_string, "%04d%02d%02d%02d%02d%02d", pnow->tm_year + 1900, pnow->tm_mon+1, pnow->tm_mday, pnow->tm_hour, pnow->tm_min, pnow->tm_sec);

    return getFileNameWithoutExt(input_file) + '_' + ((string)time_string) + ".dat";
}

void
printParts(const vector<particle> &parts, ostream &st)
{
    for(int i = 0; i < parts.size(); ++i) {
        st << scientific << setprecision(5)
             << parts[i].m << " " << parts[i].x << " " << parts[i].y << " " << parts[i].z << " "
                                  << parts[i].vx << " " << parts[i].vy << " " << parts[i].vz << endl;
    }
}

int
main(int argc, char *argv[])
{
    vector<string> args(argv, argv + argc);
    size_t N, count;
    float t, dt, init_t;
    string input_file_name = args.at(2);
    ifstream ifs(input_file_name);
    ofstream ofs(createFileName(input_file_name));

    float *buf = (float *) aligned_alloc(32, sizeof(float) * 3 * 8);

    ifs >> N >> init_t;
    if(N % 8 != 0) {
        cerr << "N is not correct!" << endl;
        return -1;
    }

    vector<particle> parts(N);

    ofs << scientific << setprecision(8);

    dt = 1000;
    for(int i = 0; i < N; ++i) {
        float m, x, y, z, vx, vy, vz;
        ifs >> m >> x >> y >> z >> vx >> vy >> vz;
        parts[i] = {x, y, z, m, vx, vy, vz, eps2};
        dt = min(dt, eps / sqrt(vx*vx + vy*vy + vz*vz));
    }

    cout << dt << endl;

    /* print initial info */
    ofs << t << endl;
    printParts(parts, ofs);

    _mm256_zeroall();
    __m256 DT = _mm256_set1_ps(dt);
    t = 0.0;
    count = 0;
    omp_set_num_threads(stoi(args.at(1)));
    #pragma omp parallel for
    for(int i = 0; i < N; ++i) {
        parts.at(i).x += parts.at(i).vx * 0.5 * dt;
        parts.at(i).y += parts.at(i).vy * 0.5 * dt;
        parts.at(i).z += parts.at(i).vz * 0.5 * dt;
    }

    while(t < MAX_T) {
#ifdef DEBUG_TIME
        auto start = chrono::system_clock::now();
#endif
        cout << "steps: " << count << endl;
        #pragma omp parallel for
        for(int i = 0; i < N; i += 4) {
            /* Hello SIMD Registers */
            __m256 XI, YI, ZI, EPS2;
            __m256 XJ, YJ, ZJ, MJ, RJ;

            __m256 lift1, lift2;
            __m256 DX, DY, DZ;
            __m256 PHI;
            __m256 AX, AY, AZ;

            XI = _mm256_set_ps(
                        parts.at(i+3).x,
                        parts.at(i+2).x,
                        parts.at(i+1).x,
                        parts.at(i+0).x,
                        parts.at(i+3).x,
                        parts.at(i+2).x,
                        parts.at(i+1).x,
                        parts.at(i+0).x
                 );
            YI = _mm256_set_ps(
                        parts.at(i+3).y,
                        parts.at(i+2).y,
                        parts.at(i+1).y,
                        parts.at(i+0).y,
                        parts.at(i+3).y,
                        parts.at(i+2).y,
                        parts.at(i+1).y,
                        parts.at(i+0).y
                 );
            ZI = _mm256_set_ps(
                        parts.at(i+3).z,
                        parts.at(i+2).z,
                        parts.at(i+1).z,
                        parts.at(i+0).z,
                        parts.at(i+3).z,
                        parts.at(i+2).z,
                        parts.at(i+1).z,
                        parts.at(i+0).z
                 );
            EPS2 = _mm256_set_ps(
                        parts.at(i+3).eps2,
                        parts.at(i+2).eps2,
                        parts.at(i+1).eps2,
                        parts.at(i+0).eps2,
                        parts.at(i+3).eps2,
                        parts.at(i+2).eps2,
                        parts.at(i+1).eps2,
                        parts.at(i+0).eps2
                  );
            AX = _mm256_setzero_ps();
            AY = _mm256_setzero_ps();
            AZ = _mm256_setzero_ps();
            for(int j = 0; j < N; j += 2) {
                lift1 = _mm256_load_ps(parts.at(j).f);
                lift2 = _mm256_load_ps(parts.at(j+1).f);

                /* XJ: MSB {{xj1, yj1, zj1, mj1}, {xj0, yj0, zj0, mj0}} */
                XJ = _mm256_insertf128_ps(lift1, _mm256_castps256_ps128(lift2), 1);

                MJ = _mm256_shuffle_ps(XJ, XJ, 0b11111111);
                YJ = _mm256_shuffle_ps(XJ, XJ, 0b01010101);
                ZJ = _mm256_shuffle_ps(XJ, XJ, 0b10101010);
                XJ = _mm256_shuffle_ps(XJ, XJ, 0b00000000);

                DX = _mm256_sub_ps(XJ, XI);
                DY = _mm256_sub_ps(YJ, YI);
                DZ = _mm256_sub_ps(ZJ, ZI);

                /* r^2 = dx^2+dy^2+dz^2+eps^2 */
                RJ = _mm256_setzero_ps();
                RJ = _mm256_fmadd_ps(DX, DX, RJ);
                RJ = _mm256_fmadd_ps(DY, DY, RJ);
                RJ = _mm256_fmadd_ps(DZ, DZ, RJ);
                RJ = _mm256_add_ps(RJ, EPS2);

                /* culculate m/r^3 */
                RJ = _mm256_rsqrt_ps(RJ);
                MJ = _mm256_mul_ps(RJ, MJ);
                RJ = _mm256_mul_ps(RJ, RJ);
                RJ = _mm256_mul_ps(RJ, MJ);

                AX = _mm256_fmadd_ps(DX, RJ, AX);
                AY = _mm256_fmadd_ps(DY, RJ, AY);
                AZ = _mm256_fmadd_ps(DZ, RJ, AZ);
            }
            AX = _mm256_mul_ps(AX, DT);
            DX = _mm256_castps128_ps256(_mm256_extractf128_ps(AX, 1));
            AX = _mm256_add_ps(AX, DX);
            _mm256_store_ps(buf+0, AX);

            AY = _mm256_mul_ps(AY, DT);
            DY = _mm256_castps128_ps256(_mm256_extractf128_ps(AY, 1));
            AY = _mm256_add_ps(AY, DY);
            _mm256_store_ps(buf+8, AY);

            AZ = _mm256_mul_ps(AZ, DT);
            DZ = _mm256_castps128_ps256(_mm256_extractf128_ps(AZ, 1));
            AZ = _mm256_add_ps(AZ, DZ);
            _mm256_store_ps(buf+16, AZ);

            for(int k = 0; k < 4; ++k) {
                parts.at(i+k).vx += buf[0+k];
                parts.at(i+k).vy += buf[8+k];
                parts.at(i+k).vz += buf[16+k];
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            parts.at(i).x += parts.at(i).vx * dt;
            parts.at(i).y += parts.at(i).vy * dt;
            parts.at(i).z += parts.at(i).vz * dt;
        }
        t += dt;
        count++;
        if(count % 1000 == 0) {
            ofs << t << endl;
            printParts(parts, ofs);
        }
#ifdef DEBUG_TIME
        auto end = chrono::system_clock::now();
        cout << "time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
#endif
    }
    cout << "total run: " << count << endl;
    return 0;
}
