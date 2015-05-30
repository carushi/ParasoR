#ifndef HOMEDIR
#define HOMEDIR "/home/"
#endif

#ifndef _PARAM_HH
#define _PARAM_HH

#include <vector>
#include <cstdlib>
#include <string>
#include "convert.hh"
#include "sequence.hh"

namespace Rfold {
namespace Parameter {

using std::cerr;

class Convert;

extern DOUBLE loghairpin[MAXLOOP+1];
// extern DOUBLE logtetra[MAXLOOP+1];
extern DOUBLE logmismatchH[7][5][5];
extern DOUBLE logmismatchI[7][5][5];
extern DOUBLE logmismatchM[7][5][5];
extern DOUBLE logmismatch1nI[7][5][5];
extern DOUBLE logmismatch23I[7][5][5];
extern DOUBLE logmismatchExt[7][5][5];
extern DOUBLE Triloop[40];
extern DOUBLE Tetraloop[40];
extern DOUBLE Hexaloop[40];
extern DOUBLE logstack[7][7];
extern DOUBLE logbulge[MAXLOOP+1];
extern DOUBLE logTermAU;
extern DOUBLE logint11[8][8][5][5];
extern DOUBLE logint21[8][8][5][5][5];
extern DOUBLE logint22[8][8][5][5][5][5];
extern DOUBLE loginternal[31];
extern DOUBLE logdangle5[8][5];
extern DOUBLE logdangle3[8][5];
extern DOUBLE logninio[MAXLOOP+1];
extern bool initialized;
extern bool inittermau;
extern DOUBLE logMLintern;
extern DOUBLE logMLclosing;
extern DOUBLE logML_BASE;
extern DOUBLE lxc37;

extern std::string Triloops;
extern std::string Tetraloops;
extern std::string Hexaloops;

static const bool no_closingGU = false;
static const bool tetra = true;

inline bool IsAU(const int type) {
    return (type > 2);
}
inline bool IsCloseGU(const int type) {
    return (type == 3 || type == 4);
}

static void SetTemperature(int temp)
{
    ChangeTemperature(temp);
    cout << "-Temperature " << temp << " " << kT << endl;
}

static void PrintSummary()
{
    cout << "Hairpin score:\t" << ExpEnergy(loghairpin[3]) << " " << ExpEnergy(loghairpin[4]) << " ...\n";
    cout << "Bulge score:\t" << ExpEnergy(logbulge[0]) << " " << ExpEnergy(logbulge[1]) << " ...\n";
    cout << "Internal score:\t" << ExpEnergy(loginternal[2]) << " " << ExpEnergy(loginternal[3]) << " ...\n";
    cout << "Tri:\t" << ExpEnergy(Triloop[0]) << " " << ExpEnergy(Triloop[1]) << " ...\n";
    cout << "Tetra:\t" << ExpEnergy(Tetraloop[0]) << " " << ExpEnergy(Tetraloop[1]) << " ...\n";
    cout << "TermAU:\t" << ExpEnergy(logTermAU) << "\n";
    cout << "Ninio:\t" << ExpEnergy(logninio[0]) << " " << ExpEnergy(logninio[1]) << " ...\n";
    cout << "ML_internal:\t" << ExpEnergy(logMLintern) << "\n";
    cout << "ML_closing:\t" << ExpEnergy(logMLclosing) << "\n";
    cout << "ML_BASE:\t" << ExpEnergy(logML_BASE) << "\n";
    cout << "Triloops:\t" << Triloops << "\n";
    cout << "Tetraloops:\t" << Tetraloops << "\n";
    cout << "Hexaloops:\t" << Hexaloops << "\n";
    cout << "Lxc:\t" << lxc37 << endl;

}
static void ChangeEnergyParam(string name = "")
{
    bool doubt = false;
    string file;
    if (name.length() == 0 || name == "Turner2004" || name[0] == 'T')
        file = HOMEDIR+string("energy_param/rna_turner2004.par");
    else if (name == "Andronescu" || name[0] == 'A')
        file = HOMEDIR+string("energy_param/rna_andronescu2007.par");
    else if (name == "Turner1999")
        file = HOMEDIR+string("energy_param/rna_turner1999.par");
    else {
        file = name;
        doubt = true;
    }
    cout << "-Read Energy " << name << " -> " << file << endl;
    class Convert convert;
    if (!convert.ConvertParamFile(file))
        cerr << "Format error: Energy param file " << file << endl;
    else
        initialized = true;
    if (doubt) PrintSummary();
}

static void InitEnergyParam()
{
    if (!initialized) ChangeEnergyParam();
}

/**
 * @return Dangling energy.
 * Based on Vienna RNA package 2.0.7.
 */
// static DOUBLE SumDangle(int type, LEN five, LEN three, const Sequence& seq)
// {
//     DOUBLE temp = 0.0;
//     if (five > 0) temp = Logsum(temp, logdangle5[type][seq.seqget(five)]);
//     if (three <= (int)seq.length) temp = Logsum(temp, logdangle3[type][seq.seqget(three)]);
//     else if (IsAU(type)) temp = Logsum(temp, logTermAU);
//     return temp;
// }

/**
 * @param five Base type dangled at 5'.
 * @param three Base type dangled at 3'.
 * @param ext use Mismatch External parameters instead of mismatch Multi.
 * @return External or normal Multi loop stem energy.
 * From exp_E_MLstem in Vienna RNA package 2.0.7.
 */

static DOUBLE SumExtML(int type, LEN five, LEN three, bool ext, const Sequence& seq)
{
    DOUBLE temp = 0.0;
    if (five > 0 && three <= seq.length) {
        if (ext) temp = Logsum(temp, logmismatchExt[type][seq.seqget(five)][seq.seqget(three)]);
        else temp = Logsum(temp, logmismatchM[type][seq.seqget(five)][seq.seqget(three)]);
        if (IsAU(type)) temp = Logsum(temp, logTermAU);
    } else {
        if (five > 0) temp = Logsum(temp, logdangle5[type][seq.seqget(five)]);
        if (three <= seq.length) temp = Logsum(temp, logdangle3[type][seq.seqget(three)]);
        if (IsAU(type)) temp = Logsum(temp, logTermAU);
    }
    return temp;
}

/**
 * @return Hairpin energy with base pair between 'i' and 'j'.
 * From exp_E_Hairpin in Vienna RNA package 2.0.7.
 */
static DOUBLE LogHairpinEnergy(LEN i, LEN j, const Sequence& seq)
{
    int type = seq.slidebp(i, j);
    DOUBLE q = -INF;
    LEN d = j-i-1;
    q = (d <= MAXLOOP) ? loghairpin[d] : loghairpin[MAXLOOP]-lxc37*log(d/(DOUBLE)MAXLOOP)*10.0/kT;
    if (d < 3) return q;
    if (tetra && d==4) {
        string sub_seq = seq.substr(i,d+2);
        size_t tel = Tetraloops.find(sub_seq);
        if (tel != string::npos) {
            if (type != 7) return Tetraloop[tel/7];
            else q = Logsum(q, Tetraloop[tel/7]);
        }
    }
    if (tetra && d==6) {
        string sub_seq = seq.substr(i,d+2);
        size_t tel = Hexaloops.find(sub_seq);
        if (tel != string::npos) return Hexaloop[tel/9];
    }
    if (d == 3) {
        string sub_seq = seq.substr(i,d+2);
        size_t tel = Triloops.find(sub_seq);
        if (tel != string::npos) return Triloop[tel/6];
        if (IsAU(type)) q = Logsum(q, logTermAU);
    } else
        q = Logsum(q, logmismatchH[type][seq.seqget(i+1)][seq.seqget(j-1)]);
    return q;
}

/**
 * @param u Length of bulge.
 * @return Bulge loop energy surrounded by base pair of 'type' and 'type2'.
 * From exp_E_IntLoop in Vienna RNA package 2.0.7.
 */
static DOUBLE LogBulge(LEN u, int type, int type2)
{
    DOUBLE z = logbulge[u];
    if (u == 1) z = Logsum(z, logstack[type][type2]);
    else {
        if (IsAU(type)) z = Logsum(z, logTermAU);
        if (IsAU(type2)) z = Logsum(z, logTermAU);
    }
    return z;
}

/**
 * Calculates loop energy (Must satisfy i <= p < q <= j).
 * From exp_E_IntLoop in Vienna RNA package 2.0.7.
 * @return Internal loop energy surrounded by base pair between 'i' and 'j', and 'p' and 'q'.
 */
static DOUBLE LogLoopEnergy(LEN i, LEN j, LEN p, LEN q, const Sequence&  seq)
{
    int type = seq.slidebp(i, j), type2 = seq.sliderbp(p, q);
    DOUBLE z = -INF;
    LEN u1 = p-i-1, u2 = j-q-1, u = max(u1, u2);
    if (u1 == 0 && u2 == 0) {
        return logstack[type][type2];
    } else if (no_closingGU && (IsCloseGU(type) || IsCloseGU(type2))) {
        return z;
    } else if ((u1 == 0) || (u2 == 0)) { /* bulge */
        return LogBulge(u, type, type2);
    } else {
        if (u <= 2) {             /* short internal */
            if (u1+u2 == 2) {
                z = logint11[type][type2][seq.seqget(i+1)][seq.seqget(j-1)];
            } else if (u1 == 1 && u2 == 2) {
                z = logint21[type][type2][seq.seqget(i+1)][seq.seqget(q+1)][seq.seqget(j-1)];
            } else if (u1 == 2 && u2 == 1) {
                z = logint21[type2][type][seq.seqget(q+1)][seq.seqget(i+1)][seq.seqget(p-1)];
            } else {
                z = logint22[type][type2][seq.seqget(i+1)][seq.seqget(p-1)][seq.seqget(q+1)][seq.seqget(j-1)];
            }
        } else {                  /* long internal */
            z = loginternal[u1+u2];
            DOUBLE temp1, temp2;
            if (u1 == 1 || u2 == 1) {
                temp1 = logmismatch1nI[type][seq.seqget(i+1)][seq.seqget(j-1)];
                temp2 = logmismatch1nI[type2][seq.seqget(q+1)][seq.seqget(p-1)];
            } else if (u1+u2 == 5) {
                temp1 = logmismatch23I[type][seq.seqget(i+1)][seq.seqget(j-1)];
                temp2 = logmismatch23I[type2][seq.seqget(q+1)][seq.seqget(p-1)];
            } else {
                temp1 = logmismatchI[type][seq.seqget(i+1)][seq.seqget(j-1)];
                temp2 = logmismatchI[type2][seq.seqget(q+1)][seq.seqget(p-1)];
            }
            DOUBLE temp3 = logninio[abs(u1-u2)];
            z = Logsum(z, temp1, temp2, temp3);
        }
    }
    return z;
}
}
}

#endif