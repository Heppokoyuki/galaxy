#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <omp.h>
#include <sstream>

const double eps = 0.01;
const double MAX_T = 10;

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
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

template<typename T>
inline
T
square(T e) {
    return e * e;
}

double
calcDistance2(const double v1[3], const double v2[3])
{
    return square(v1[0] - v2[0]) + square(v1[1] - v2[1]) + square(v1[2] - v2[2]);
}

double
calcPotential(const vector<particle> &parts, int i, int b, const vector<vector<double>> &dist)
{
    double res = 0.0;
    for(int j = 0; j < parts.size(); ++j) {
        res += parts[j].m * (parts[j].x[b] - parts[i].x[b]) / pow(dist[j][i] + eps*eps, 1.5);
    }
    return res;
}

void
printParts(const vector<particle> &parts, ostream &st)
{
    for(int i = 0; i < parts.size(); ++i) {
        st << scientific << setprecision(5)
             << parts[i].m << " " << parts[i].x[0] << " " << parts[i].x[1] << " " << parts[i].x[2] << " "
                                  << parts[i].v[0] << " " << parts[i].v[1] << " " << parts[i].v[2] << endl;
    }
}

int main()
{
    size_t N, count;
    double t, dt, init_t;
    string input_file_name = "unisp4k.dat";
    ifstream ifs(input_file_name);
    ofstream ofs(createFileName(input_file_name));

    ifs >> N >> init_t;

    vector<particle> parts(N);
    vector<vector<double>> distance(N, vector<double>(N));

    ofs << scientific << setprecision(8);

    dt = 1000;
    for(int i = 0; i < N; ++i) {
        double m, x, y, z, vx, vy, vz;
        ifs >> m >> x >> y >> z >> vx >> vy >> vz;
        parts[i] = {m, {x, y, z}, {vx, vy, vz}};
        dt = min({dt, abs(eps/vx), abs(eps/vy), abs(eps/vz)});
    }

    cout << dt << endl;

    t = 0.0;
    count = 0;
    /* print initial info */
    ofs << t << endl;
    printParts(parts, ofs);

    for(int idx = 0; idx < N; ++idx) {
        for(int base = 0; base < 3; ++base) {
            parts[idx].x[base] = parts[idx].x[base] + parts[idx].v[base] * 0.5 * dt;
        }
    }

    while(t < MAX_T) {
        cout << "steps: " << count << endl;

        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                distance[i][j] = calcDistance2(parts[i].x, parts[j].x);
            }
        }

        #pragma omp parallel for
        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].v[base] = parts[idx].v[base] + calcPotential(parts, idx, base, distance) * dt;
            }
        }

        #pragma omp parallel for
        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].x[base] = parts[idx].x[base] + parts[idx].v[base] * dt;
            }
        }

        t += dt;
        count++;
        if(count % 10 == 0) {
            ofs << N << t << endl;
            printParts(parts, ofs);
        }
    }

    cout << "total run: " << count << endl;
    return 0;
}
