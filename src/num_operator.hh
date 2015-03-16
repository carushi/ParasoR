#ifndef _NUM_OPERATOR_HH
#define _NUM_OPERATOR_HH

#include "energy_const.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <climits>

#define SCALE (10.0)

#define SMOOTH(X) ((X)/SCALE < -1.2283697) ? 0.0 : (((X)/SCALE > 0.8660254) ? (X) : SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1.0)*(sin((X)/SCALE-0.34242663)+1.0))
//#define SMOOTH(X) ((X) < 0 ? 0 : (X))

extern DOUBLE kT;
extern const int temperature;

inline DOUBLE Correlation(const std::vector<DOUBLE>& ori, const std::vector<DOUBLE>& mut)
{
    DOUBLE sum_sq_x = 0.0, sum_sq_y = 0.0, sum_coproduct = 0.0, mean_x = 0.0, mean_y = 0.0, n = 1.0;
    for (std::vector<DOUBLE>::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++)
    {
        DOUBLE x = *it, y = *it2;
        DOUBLE delta_x = (x-mean_x), delta_y = (y-mean_y), sweep = (n-1.0)/n;
        sum_sq_x += delta_x * delta_x * sweep;
        sum_sq_y += delta_y * delta_y * sweep;
        sum_coproduct += delta_x * delta_y * sweep;
        mean_x += delta_x / n;
        mean_y += delta_y / n++;
    }
    if (sum_sq_x == 0.0 || sum_sq_y == 0.0) return std::numeric_limits<DOUBLE>::quiet_NaN();
    else return sum_coproduct/sqrt(sum_sq_x*sum_sq_y);
}

inline static bool Is_INF(const DOUBLE value) {
    return (value <= -INF);
}

inline static DOUBLE Logsumexp(DOUBLE x, DOUBLE y) {
    if (x <= -INF) return y;
    else if (y <= -INF) return x;
    else if (x > y) return (x + log1p(exp(-x+y)));
    else return (y + log1p(exp(-y+x)));
}

inline static DOUBLE Logsum(DOUBLE a, DOUBLE b) {
    if (Is_INF(a) || Is_INF(b)) return -INF;
    else return a+b;
}
inline static DOUBLE Logsum(DOUBLE a, DOUBLE b, DOUBLE c) {
    return Logsum(a, Logsum(b, c));
}
inline static DOUBLE Logsum(DOUBLE a, DOUBLE b, DOUBLE c, DOUBLE d) {
    return Logsum(Logsum(a, b), Logsum(c, d));
}

static void ChangeTemperature(const int value) {
    kT = (value+K0)*GASCONST;
}
static DOUBLE LogEnergy(const int value) {
    if (value == INTINF) return -INF;
    else return -(static_cast<DOUBLE>(value*10.0))/kT;
}
static DOUBLE LogSmoothingEnergy(const int value) {
    if (value == INTINF) return -INF;
    else return SMOOTH(-(DOUBLE)value)*10.0/kT;
}
static DOUBLE ExpEnergy(const DOUBLE value) {
    if (Is_INF(value)) return INTINF;
    else return static_cast<DOUBLE>(-value/10.0*kT);
}
static DOUBLE Energy(const DOUBLE value) {
    return -kT*(value)/1000.;
}
#endif


