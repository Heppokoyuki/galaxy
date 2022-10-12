#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

const double eps = 1e-5;
const double MAX_T = 1e-4;

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
};

int main()
{
    int N, count;
    double t, dt, init_t;
    ifstream ifs("pl8k.dat");
    ofstream ofs("output.dat");
    ifs >> N >> init_t;
    vector<particle> parts(N);

    ofs << scientific << setprecision(5);

    /* calc the potential energy */
    auto potential = [N, &parts](int i, int b) {
        double res = 0.0;
        for(int j = 0; j < N; ++j) {
            res += parts[j].m * (parts[j].x[b] - parts[i].x[b]) / pow(pow((parts[j].x[b] - parts[i].x[b]), 2) + eps*eps, 3/2);
        }
        return res;
    };

    auto printParts = [N, &parts](ostream &st) {
        for(int i = 0; i < N; ++i) {
            st << scientific << setprecision(5)
                 << parts[i].m << " " << parts[i].x[0] << " " << parts[i].x[1] << " " << parts[i].x[2] << " "
                                      << parts[i].v[0] << " " << parts[i].v[1] << " " << parts[i].v[2] << endl;
        }
    };

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
        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].x[base] = parts[idx].x[base] + parts[idx].v[base] * dt;
            }
        }

        for(int idx = 0; idx < N; ++idx) {
            for(int base = 0; base < 3; ++base) {
                parts[idx].v[base] = parts[idx].v[base] + potential(idx, base) * dt;
            }
        }
        printParts(ofs);

        t += dt;
        count++;
    }

    printParts(cout);
    cout << "total run: " << count << endl;
    return 0;
}
