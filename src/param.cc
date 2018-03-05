
#include "param.hh"

namespace Rfold {
namespace Parameter {

/**
 * Parameters. Based on Vienna RNA package 2.0.7.
 */
DOUBLE loghairpin[MAXLOOP+1];
// DOUBLE logtetra[MAXLOOP+1];
DOUBLE logmismatchH[NBPAIRS+1][5][5];
DOUBLE logmismatchI[NBPAIRS+1][5][5];
DOUBLE logmismatchM[NBPAIRS+1][5][5];
DOUBLE logmismatch1nI[NBPAIRS+1][5][5];
DOUBLE logmismatch23I[NBPAIRS+1][5][5];
DOUBLE logmismatchExt[NBPAIRS+1][5][5];
DOUBLE Triloop[40];
DOUBLE Tetraloop[40];
DOUBLE Hexaloop[40];
DOUBLE logstack[NBPAIRS+1][NBPAIRS+1];
DOUBLE logbulge[MAXLOOP+1];

DOUBLE logint11[NBPAIRS+1][NBPAIRS+1][5][5];
DOUBLE logint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
DOUBLE logint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
DOUBLE loginternal[31];
DOUBLE logdangle5[NBPAIRS+1][5];
DOUBLE logdangle3[NBPAIRS+1][5];
DOUBLE logninio[MAXLOOP+1];
DOUBLE logMLintern;
DOUBLE logMLclosing;
DOUBLE logML_BASE;
DOUBLE logTermAU;

std::string Triloops;
std::string Tetraloops;
std::string Hexaloops;

DOUBLE lxc37 = 107.856;
bool initialized = false;
bool inittermau = false;
bool old_param = false;


}
}
