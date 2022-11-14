#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>

const size_t N = 4096;
const double dr = 1e-2;
const double eps = 0.01;

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
};

double
calcx2(const particle &tmp) {
    return pow(tmp.x[0], 2) + pow(tmp.x[1], 2) + pow(tmp.x[2], 2);
}

double
calcv2(const particle &tmp) {
    return pow(tmp.v[0], 2) + pow(tmp.v[1], 2) + pow(tmp.v[2], 2);
}

double
distance(const double v1[3], const double v2[3])
{
    return sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2) + eps*eps);
}

double
calcPotential(const vector<particle> &parts, const int steps, int i)
{
    double res = 0.0;
    i += N * steps;
    for(int j = N * steps; j < N * (steps+1); ++j) {
        if(i == j) continue;
        res += parts[j].m / distance(parts[j].x, parts[i].x);
    }
    return -1.0 * res;
}

void
printMassDist(const double &t, const vector<particle> &v, string filename)
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
        r.push_back(calcx2(v[i]));
    }

    moment = 0;
    som = 0;
    for(int i = 0; i < N; ++i) {
        tmp = v[i];
        moment += tmp.m * calcx2(tmp);
        som += tmp.m;
    }
    cog = moment / som;
    for(int i = 0; i < N; ++i) {
        r[i] = abs(r[i] - cog);
    }

    sort(r.begin(), r.end());

    output << scientific << setprecision(8);
    output << t << endl;
    index = 0;
    for(int i = 0; i < 10000; ++i) {
        for(; index < N; ++index) {
            if(r[index] > (double)i * dr) break;
        }
        output << i * dr << " " << index << endl;
    }
}

void
mass(const vector<double> &t, const vector<particle> &v) {
    vector<particle> m[10];
    size_t steps = t.size();
    int factor = (int)steps / 10;

    for(int i = 0; i < 10; ++i) {
        for(int j = 0; j < N; ++j) {
            m[i].push_back(v[factor * i * N + j]);
        }
        printMassDist(t[factor * i], m[i], "mass" + to_string(i) + ".dat");
    }
}

void
energy(const vector<double> &t, const vector<particle> &v) {
    ofstream output1("penergy.dat"), output2("kenergy.dat");
    size_t steps;
    particle tmp;
    double kenergy, penergy;

    steps = t.size() - 1;
    output1 << scientific << setprecision(8);
    output2 << scientific << setprecision(8);

    for(int i = 0; i < steps; ++i) {
        kenergy = 0;
        penergy = 0;

        for(int j = 0; j < N; ++j) {
            tmp = v[i * N + j];
            kenergy += 0.5 * tmp.m * calcv2(tmp);
            penergy += 0.5 * tmp.m * calcPotential(v, i, j);
        }
        output1 << t[i] << " " << penergy << endl;
        output2 << t[i] << " " << kenergy << endl;
    }
}


int main()
{
    ifstream ifs("output.dat");
    vector<particle> parts;
    vector<double> time;
    double t;
    particle tmp;

    do {
        ifs >> t;
        time.push_back(t);
        for(int i = 0; i < N; ++i) {
            ifs >> tmp.m >> tmp.x[0] >> tmp.x[1] >> tmp.x[2] >> tmp.v[0] >> tmp.v[1] >> tmp.v[2];
            parts.push_back(tmp);
        }
    }
    while(!ifs.eof());

    energy(time, parts);
    mass(time, parts);

    return 0;
}
