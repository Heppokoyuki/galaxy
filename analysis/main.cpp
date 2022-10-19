#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>

const size_t N = 4096;
const double dr = 1e-2;
const double eps = 1e-5;

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
};

double
calcx(const particle &tmp) {
    return pow(tmp.x[0], 2) + pow(tmp.x[1], 2) + pow(tmp.x[2], 2);
}

double
calcv2(const particle &tmp) {
    return pow(tmp.v[0], 2) + pow(tmp.v[1], 2) + pow(tmp.v[2], 2);
}

double
calcPotential(const vector<particle> &parts, int steps, int i)
{
    double res = 0.0;
    i += N * steps;
    for(int j = N * steps; j < N * (steps+1); ++j) {
        for(int b = 0; b < 3; ++b) {
            res += parts[j].m * (parts[j].x[b] - parts[i].x[b]) / pow(pow((parts[j].x[b] - parts[i].x[b]), 2) + eps*eps, 3/2);
        }
    }
    return res;
}

void
printMassDist(const vector<particle> &v, string filename)
{
    ofstream output(filename);
    particle tmp;
    vector<double> r;
    int index;
    double moment, som, cog;

    if(v.size() != N) {
        cout << "vector size is not correct!" << endl;
    }

    for(int i = 0; i < N; ++i) {
        r.push_back(calcx(v[i]));
    }

    moment = 0;
    som = 0;
    for(int i = 0; i < N; ++i) {
        tmp = v[i];
        moment += tmp.m * calcx(tmp);
        som += tmp.m;
    }
    cog = moment / som;
    for(int i = 0; i < N; ++i) {
        r[i] = abs(r[i] - cog);
    }

    sort(r.begin(), r.end());

    output << scientific << setprecision(8);
    index = 0;
    for(int i = 0; i < 10000; ++i) {
        for(; index < N; ++index) {
            if(r[index] > (double)i * dr) break;
        }
        output << i * dr << " " << index << endl;
    }
}

void
mass(const vector<particle> &v) {
    vector<particle> m[10];

    for(int i = 0; i < 10; ++i) {
        for(int j = 0; j < N; ++j) {
            m[i].push_back(v[i * 400 * N + j]);
        }
        printMassDist(m[i], "mass" + to_string(i) + ".dat");
    }
}

void
energy(vector<particle> &v) {
    ofstream output("energy.dat");
    size_t steps;
    double t = 0;
    double dt = 6.75493e-06;
    particle tmp;
    double kenergy, penergy;

    steps = v.size() / N;
    output << scientific << setprecision(8);

    for(int i = 0; i < steps; ++i) {
        kenergy = 0;
        penergy = 0;
        t += dt;
        for(int j = 0; j < N; ++j) {
            tmp = v[i * N + j];
            kenergy += 0.5 * tmp.m * calcv2(tmp);
            penergy += calcPotential(v, i, j);
        }
        output << t << " " << kenergy + penergy * 0.5 << endl;
    }
}


int main()
{
    ifstream ifs("4701.dat");
    vector<particle> parts;
    particle tmp;

    while(!ifs.eof()) {
        ifs >> tmp.m >> tmp.x[0] >> tmp.x[1] >> tmp.x[2] >> tmp.v[0] >> tmp.v[1] >> tmp.v[2];
        parts.push_back(tmp);
    }

    energy(parts);
    mass(parts);

    return 0;
}
