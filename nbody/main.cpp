#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <omp.h>

const double eps = 1e-5;
const double MAX_T = 1e-3;

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
};

double
calcPotential(const vector<particle> &parts, int i, int b)
{
    double res = 0.0;
    for(int j = 0; j < parts.size(); ++j) {
        res += parts[j].m * (parts[j].x[b] - parts[i].x[b]) / pow(pow((parts[j].x[b] - parts[i].x[b]), 2) + eps*eps, 3/2);
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
    int N, count;
    double t, dt, init_t;
    ifstream ifs("pl4k.dat");
    ofstream ofs("output.dat");
    ifs >> N >> init_t;
    vector<particle> parts(N);
    vector<particle> bparts(N);

    ofs << scientific << setprecision(5);

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
    while(t < MAX_T) {
        cout << "steps: " << count << endl;
        #pragma omp parallel for
        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].x[base] = parts[idx].x[base] + parts[idx].v[base] * dt;
            }
        }

        #pragma omp parallel for
        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].v[base] = parts[idx].v[base] + calcPotential(parts, idx, base) * dt;
            }
        }

        printParts(parts, ofs);

        t += dt;
        count++;
    }

    cout << "total run: " << count << endl;
    return 0;
}
