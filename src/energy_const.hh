#ifndef _ENERGY_CONST_HH
#define _ENERGY_CONST_HH
#include <limits>
#if defined (LONG)
#define DOUBLE long double
#define INF (std::numeric_limits<long double>::max()/100.0)
#elif defined (SHORT)
#define DOUBLE float
#define INF (std::numeric_limits<float>::max()/10.0)
#else
#define DOUBLE double
#define INF (std::numeric_limits<double>::max()/100.0)
#endif

#define LEN long int

#define GASCONST (1.98717)  /* in [cal/K] */
#define K0  (273.15)
#define INTINF (1000000)
#define TURN (3)
#define MAXLOOP (30)



#endif
