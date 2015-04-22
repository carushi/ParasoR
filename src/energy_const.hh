#ifndef _ENERGY_CONST_HH
#define _ENERGY_CONST_HH

#if defined (LONG)
#define DOUBLE long double
#define INF (100000000000)
#elif defined (SHORT)
#define DOUBLE float
#define INF (100000)
#else
#define DOUBLE double
#define INF (100000000000)
#endif

#define LEN long int

#define GASCONST (1.98717)  /* in [cal/K] */
#define K0  (273.15)
#define INTINF (1000000)
#define TURN (3)
#define MAXLOOP (30)



#endif