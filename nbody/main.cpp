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

// #define DEBUG_TIME

const double eps = 0.01;
const double eps2 = eps*eps;
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

void
printParts(const vector<particle> &parts, ostream &st)
{
    for(int i = 0; i < parts.size(); ++i) {
        st << scientific << setprecision(5)
             << parts[i].m << " " << parts[i].x[0] << " " << parts[i].x[1] << " " << parts[i].x[2] << " "
                                  << parts[i].v[0] << " " << parts[i].v[1] << " " << parts[i].v[2] << endl;
    }
}

int main(int argc, char *argv[])
{
    vector<string> args(argv, argv + argc);
    size_t N, count;
    double t, dt, init_t;
    string input_file_name = args.at(2);
    ifstream ifs(input_file_name);
    ofstream ofs(createFileName(input_file_name));

    ifs >> N >> init_t;

    vector<particle> parts(N);
    double distance;

    ofs << scientific << setprecision(8);

    dt = 1000;
    for(int i = 0; i < N; ++i) {
        double m, x, y, z, vx, vy, vz;
        ifs >> m >> x >> y >> z >> vx >> vy >> vz;
        parts[i] = {m, {x, y, z}, {vx, vy, vz}};
        dt = min(dt, eps / sqrt(vx*vx + vy*vy + vz*vz));
    }

    cout << dt << endl;

    t = 0.0;
    count = 0;
    /* print initial info */
    ofs << t << endl;
    printParts(parts, ofs);
    
    omp_set_num_threads(stoi(args.at(1)));
    #pragma omp parallel for
    for(int idx = 0; idx < N; ++idx) {
        for(int base = 0; base < 3; ++base) {
            parts[idx].x[base] = parts[idx].x[base] + parts[idx].v[base] * 0.5 * dt;
        }
    }

    while(t < MAX_T) {
#ifdef DEBUG_TIME
        auto start = chrono::system_clock::now();
#endif
        cout << "steps: " << count << endl;

        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                double dx = parts[j].x[0] - parts[i].x[0];
                double dy = parts[j].x[1] - parts[i].x[1];
                double dz = parts[j].x[2] - parts[i].x[2];
                double rsq = dx*dx + dy*dy + dz*dz;
                double mrinv3 = parts[j].m * pow(rsq + eps2, -1.5) * dt;

                parts[i].v[0] += dx * mrinv3;
                parts[i].v[1] += dy * mrinv3;
                parts[i].v[2] += dz * mrinv3;
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
