
#include "param.hh"

namespace Rfold {
namespace Parameter {

/**
 * Parameters. Based on Vienna RNA package 2.0.7.
 */
DOUBLE loghairpin[MAXLOOP+1];
// DOUBLE logtetra[MAXLOOP+1];
DOUBLE logmismatchH[7][5][5];
DOUBLE logmismatchI[7][5][5];
DOUBLE logmismatchM[7][5][5];
DOUBLE logmismatch1nI[7][5][5];
DOUBLE logmismatch23I[7][5][5];
DOUBLE logmismatchExt[7][5][5];
DOUBLE Triloop[40];
DOUBLE Tetraloop[40];
DOUBLE Hexaloop[40];
DOUBLE logstack[7][7];
DOUBLE logbulge[MAXLOOP+1];

DOUBLE logint11[8][8][5][5];
DOUBLE logint21[8][8][5][5][5];
DOUBLE logint22[8][8][5][5][5][5];
DOUBLE loginternal[31];
DOUBLE logdangle5[8][5];
DOUBLE logdangle3[8][5];
DOUBLE logninio[MAXLOOP+1];
DOUBLE logMLintern;
DOUBLE logMLclosing;
DOUBLE logML_BASE;
DOUBLE logTermAU;

std::string Triloops;
std::string Tetraloops;
std::string Hexaloops;

DOUBLE lxc37 = 107.856;
bool initialized;
bool inittermau;


}
}