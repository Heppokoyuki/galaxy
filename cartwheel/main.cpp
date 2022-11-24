#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

using namespace std;

struct particle {
    double m;
    double x[3];
    double v[3];
};

class galaxy {
public:
    vector<particle> parts;
    int N;
    double init_t;
    double R;
    galaxy(string f);
    void setInitial(double x, double v);

private:
    double calcDistanceFromCentre(size_t i);
    void calcRadius();
};

template<typename T>
inline
T
sq(T e) {
    return e * e;
}

galaxy::galaxy(string f)
{
    ifstream ifs(f);
    ifs >> N >> init_t;
    parts.resize(N);
    for(int i = 0; i < N; ++i) {
        double m, x, y, z, vx, vy, vz;
        ifs >> m >> x >> y >> z >> vx >> vy >> vz;
        parts[i] = {m, {x, y, z}, {vx, vy, vz}};
    }
    calcRadius();
}

void
galaxy::setInitial(double x, double v)
{
    for(int i = 0; i < N; ++i) {
        parts[i].x[2] += x;
        parts[i].v[2] += v;
    }
}

double
galaxy::calcDistanceFromCentre(size_t i)
{
    double x, y, z;
    x = parts[i].x[0];
    y = parts[i].x[1];
    z = parts[i].x[2];

    return sqrt(sq(x) + sq(y) + sq(z));
}

void
galaxy::calcRadius()
{
    double res = 0;
    for(int i = 0; i < N; ++i) {
        res = max(res, calcDistanceFromCentre(i));
    }
    R = res;
}

int main()
{
    galaxy m31("m31.dat");
    galaxy plummer("pl4k.dat");
    ofstream output("cartwheel.dat");

    m31.setInitial(m31.R + plummer.R + 10, 30);

    output << scientific << setprecision(8);
    output << m31.N + plummer.N << " " << 0.0 << endl;
    for(int i = 0; i < m31.N; ++i) {
        output << right << m31.parts[i].m << " " <<
               m31.parts[i].x[0] << " " << m31.parts[i].x[1] << " " << m31.parts[i].x[2] << " " <<
               m31.parts[i].v[0] << " " << m31.parts[i].v[1] << " " << m31.parts[i].v[2] << endl;
    }
    for(int i = 0; i < plummer.N; ++i) {
        output << right << plummer.parts[i].m << " " <<
               plummer.parts[i].x[0] << " " << plummer.parts[i].x[1] << " " << plummer.parts[i].x[2] << " " <<
               plummer.parts[i].v[0] << " " << plummer.parts[i].v[1] << " " << plummer.parts[i].v[2] << endl;
    }

    return 0;
}
