#ifndef _NUM_OPERATOR_HH
#define _NUM_OPERATOR_HH

#include "energy_const.hh"
#include <cmath>
#include <iostream>

#define SCALE (10.0)

#define SMOOTH(X) ((X)/SCALE < -1.2283697) ? 0.0 : (((X)/SCALE > 0.8660254) ? (X) : SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1.0)*(sin((X)/SCALE-0.34242663)+1.0))
//#define SMOOTH(X) ((X) < 0 ? 0 : (X))

extern DOUBLE kT;
extern const int temperature;

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


