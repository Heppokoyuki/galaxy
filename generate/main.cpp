#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <limits>
#include <iomanip>

using namespace std;

const double PI = 4.0 * atan(1);
const int N = 4096;
const double eps = 0.01;

#define REP(a, b) for(int a = 0; a < b; ++a)

struct particle {
    double m;
    double x[3];
    double v[3];
};

template<typename T>
inline T sq(T a)
{
    return a * a;
}

double
calcX(double r, double phi, double mu)
{
    return r * cos(phi) * sqrt(1-sq(mu));
}

double
calcY(double r, double phi, double mu)
{
    return r * sin(phi) * sqrt(1-sq(mu));
}

double
calcZ(double r, double phi, double mu)
{
    return r * mu;
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
calcPotential(const vector<particle> &parts, int i)
{
    double res = 0.0;
    for(int j = 0; j < N; ++j) {
        if(i == j) continue;
        res += parts[j].m / distance(parts[j].x, parts[i].x);
    }
    return -1.0 * res;
}

int main()
{
    const double R = 1;
    const double m = 1 / (double)N;
    const double deviation_v = 1.0 / sqrt(10);

    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<> randR(0, R);
    uniform_real_distribution<> randM(-1, 1);
    uniform_real_distribution<> randP(0, 2*PI);
    normal_distribution<> randV(0.0, deviation_v);

    vector<particle> parts(N);
    ofstream ofs("unisp4k.dat");
    double r, phi, mu;
    double pene = 0;
    double kene = 0;
    double tmpv;

    ofs << scientific << setprecision(8);

    ofs << N << " " << 0.0 << endl;
    for(int i = 0; i < N; ++i) {
        r = cbrt(randR(mt));
        mu = randM(mt);
        phi = randP(mt);
        parts[i].m = m;
        parts[i].x[0] = calcX(r, phi, mu);
        parts[i].x[1] = calcY(r, phi, mu);
        parts[i].x[2] = calcZ(r, phi, mu);
        for(int j = 0; j < 3; ++j) {
            parts[i].v[j] = randV(mt);
        }
        ofs << right << parts[i].m << " " <<
               parts[i].x[0] << " " << parts[i].x[1] << " " << parts[i].x[2] << " " <<
               parts[i].v[0] << " " << parts[i].v[1] << " " << parts[i].v[2] << endl;
    }
    for(int i = 0; i < N; ++i) {
        particle tmp = parts[i];
        kene += 0.5 * tmp.m * calcv2(tmp);
        pene += 0.5 * tmp.m * calcPotential(parts, i);
    }
    cout << "K: " << kene << " W: " << pene << endl;
    cout << "2K/|W|: " << 2 * kene / abs(pene) << endl;
}
