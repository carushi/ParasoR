#include "convert.hh"

namespace Rfold {
namespace Parameter {

extern DOUBLE loghairpin[MAXLOOP+1];
// extern DOUBLE logtetra[MAXLOOP+1];
extern DOUBLE logmismatchH[NBPAIRS+1][5][5];
extern DOUBLE logmismatchI[NBPAIRS+1][5][5];
extern DOUBLE logmismatchM[NBPAIRS+1][5][5];
extern DOUBLE logmismatch1nI[NBPAIRS+1][5][5];
extern DOUBLE logmismatch23I[NBPAIRS+1][5][5];
extern DOUBLE logmismatchExt[NBPAIRS+1][5][5];
extern DOUBLE Triloop[40];
extern DOUBLE Tetraloop[40];
extern DOUBLE Hexaloop[40];
extern DOUBLE logstack[NBPAIRS+1][NBPAIRS+1];
extern DOUBLE logbulge[MAXLOOP+1];
extern DOUBLE logTermAU;
extern DOUBLE logint11[NBPAIRS+1][NBPAIRS+1][5][5];
extern DOUBLE logint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
extern DOUBLE logint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
extern DOUBLE loginternal[MAXLOOP+1];
extern DOUBLE logdangle5[NBPAIRS+1][5];
extern DOUBLE logdangle3[NBPAIRS+1][5];
extern DOUBLE logninio[MAXLOOP+1];
extern bool initialized;
extern bool inittermau;
extern DOUBLE logMLintern;
extern DOUBLE logMLclosing;
extern DOUBLE logML_BASE;

extern std::string Triloops;
extern std::string Tetraloops;
extern std::string Hexaloops;

extern const bool no_closingGU;
extern const bool tetra;
extern bool old_param;
extern DOUBLE lxc37;

void Convert::GetWords(string& str, vector<string>& words)
{
    words.clear();
    istringstream iss(str);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(words));
}

int Convert::param_type(string& str)
{
    if (str.length() == 0 || str[0] != '#') return -1;
    vector<string> words;
    GetWords(str, words);
    if (words.size() <= 1) return -1;
    string id = words[1];
    if (id == "stack" || id == "stack_energies") return Stack;
    else if (id == "hairpin") return Hairpin;
    else if (id == "bulge") return Bulge;
    else if (id == "interior" || id == "internal_loop") return Interior;
    else if (id == "mismatch_exterior") return MisE;
    else if (id == "mismatch_hairpin") return MisH;
    else if (id == "mismatch_interior") return MisI;
    else if (id == "mismatch_interior_1n") return Mis1n;
    else if (id == "mismatch_interior_23") return MisI23;
    else if (id == "mismatch_multi") return MisM;
    else if (id == "int11" || id == "int11_energies") return Int11;
    else if (id == "int21" || id == "int21_energies") return Int21;
    else if (id == "int22" || id == "int22_energies") return Int22;
    else if (id == "dangle5") return Dan5;
    else if (id == "dangle3") return Dan3;
    else if (id == "ML_params") return ML;
    else if (id == "NINIO") return Ninio;
    else if (id == "Triloops") return Tri;
    else if (id == "Tetraloops") return Tetra;
    else if (id == "Hexaloops") return Hexa;
    else if (id == "Misc") return Misc;
    return -1;
}

void Convert::GetArray(DOUBLE *array, int size, bool smooth = false)
{
    vector<string> words;
    for (int i = 0; i < size; ) {
        string str;
        if (!getline(ifs, str)) break;
        else if (str.length() < 2) break;
        GetWords(str, words);
        int prev = i;
        for (; i < size && i-prev < (int)words.size(); i++) {
            if (words[i-prev].find("/*") != string::npos) break;
            else {
                int temp;
                if (words[i-prev] == "INF") temp = INTINF;
                else if (words[i-prev] == "DEF") temp = -50;
                else temp = atoi(words[i-prev].c_str());
                array[i] = (smooth) ? LogSmoothingEnergy(temp) : LogEnergy(temp);
            }
        }
    }
}

void Convert::Read1Dim(DOUBLE *array, int dim, int shift, int post = 0) {
    GetArray(array+shift, dim-shift-post);
}

void Convert::Read2Dim(DOUBLE* array, int dim1, int dim2, int shift1, int shift2,
                                      int post1 = 0, int post2 = 0)
{
    if (shift1+shift2 == 0 && post1+post2 == 0)
        Read1Dim(array, dim1*dim2, 0);
    else {
        for (int i = shift1; i < dim1-post1; i++)
            Read1Dim(array+(i*dim2), dim2, shift2, post2);
    }
}

void  Convert::Read3Dim(DOUBLE *array, int dim1, int dim2, int dim3, int shift1, int shift2, int shift3,
                                       int post1 = 0, int post2 = 0, int post3 = 0)
{
    if (shift1+shift2+shift3 == 0 && post1+post2+post3 == 0)
        Read1Dim(array, dim1*dim2*dim3, 0);
    else {
        for (int i = shift1; i < dim1-post1; i++) {
            Read2Dim(array+(i*dim2*dim3), dim2, dim3, shift2, shift3, post2, post3);
        }
    }
}

void  Convert::Read4Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4,
                                       int shift1, int shift2, int shift3, int shift4,
                                       int post1 = 0, int post2 = 0, int post3 = 0, int post4 = 0)
{
    if (shift1+shift2+shift3+shift4 == 0 && post1+post2+post3+post4 == 0)
        Read1Dim(array, dim1*dim2*dim3*dim4, 0);
    else {
        for (int i = shift1; i < dim1-post1; i++) {
            Read3Dim(array+(i*dim2*dim3*dim4), dim2, dim3, dim4, shift2, shift3, shift4,
                                               post2, post3, post4);
        }
    }

}

void  Convert::Read5Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4, int dim5,
                                       int shift1, int shift2, int shift3, int shift4, int shift5,
                                       int post1 = 0, int post2 = 0, int post3 = 0, int post4 = 0, int post5 = 0)
{
    if (shift1+shift2+shift3+shift4+shift5 == 0 && post1+post2+post3+post4+post5 == 0)
        Read1Dim(array, dim1*dim2*dim3*dim4*dim5, 0);
    else {
        for (int i = shift1; i < dim1-post1; i++) {
            Read4Dim(array+(i*dim2*dim3*dim4*dim5), dim2, dim3, dim4, dim5,
                                                    shift2, shift3, shift4, shift5,
                                                    post2, post3, post4, post5);
        }
    }
}

void Convert::Read6Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                                      int shift1, int shift2, int shift3, int shift4, int shift5, int shift6,
                                      int post1 = 0, int post2 = 0, int post3 = 0, int post4 = 0, int post5 = 0, int post6 = 0)
{
    if (shift1+shift2+shift3+shift4+shift5+shift6 == 0 && post1+post2+post3+post4+post5+post6 == 0)
        Read1Dim(array, dim1*dim2*dim3*dim4*dim5*dim6, 0);
    else {
        for (int i = shift1; i < dim1-post1; i++) {
            Read5Dim(array+(i*dim2*dim3*dim4*dim5*dim6), dim2, dim3, dim4, dim5, dim6,
                                                         shift2, shift3, shift4, shift5, shift6,
                                                         post2, post3, post4, post5, post6);
        }
    }
}

void Convert::Read2DimSmooth(DOUBLE *array, int dim1, int dim2, int shift1, int shift2, int post1 = 0, int post2 = 0)
{
    if (shift1+shift2 == 0 && post1+post2 == 0)
        GetArray(array, dim1*dim2, true);
    else {
        for (int i = shift1; i < dim1-post1; i++)
            GetArray(array+(i*dim2)+shift2, dim2-shift2-post2, true);
    }
}

void Convert::Read3DimSmooth(DOUBLE *array, int dim1, int dim2, int dim3, int shift1, int shift2, int shift3,
                                       int post1 = 0, int post2 = 0, int post3 = 0)
{
    if (shift1+shift2+shift3 == 0 && post1+post2+post3 == 0)
        GetArray(array, dim1*dim2*dim3, true);
    else {
        for (int i = shift1; i < dim1-post1; i++) {
            Read2DimSmooth(array+(i*dim2*dim3), dim2, dim3, shift2, shift3, post2, post3);
        }
    }
}

void Convert::Readninio()
{
    string str;
    vector<string> words, pwords;
    int F_ninio37, MAX_NINIO;
    while(getline(ifs, str)) {
        if (str == "") break;
        GetWords(str, words);
        if (str.find("*") != string::npos) {
            pwords = words;
            continue;
        }
        if ((int)pwords.size() == (int)words.size()+2) {
            assert((int)words.size() >= 2);
            for (int i = 0; i < (int)words.size(); i++) {
                int value = atoi(words[i].c_str());
                if (pwords[i+1] == "m") F_ninio37 = value;
                else if (pwords[i+1] == "max") MAX_NINIO = value;
            }
        } else {
            assert((int)words.size() > 2);
            F_ninio37 = atoi(words[0].c_str());
            MAX_NINIO = atoi(words[2].c_str());
        }
        break;
    }
    for (int i = 0; i <= MAXLOOP; i++)
        logninio[i]=LogEnergy(min(MAX_NINIO,i*F_ninio37));
}

void Convert::ReadML()
{
    string str;
    vector<string> words, pwords;
    while(getline(ifs, str)) {
        if (str == "") break;
        GetWords(str, words);
        if (str.find("*") != string::npos) {
            pwords = words;
            continue;
        }
        if ((int)pwords.size() == (int)words.size()+2) {
            assert((int)words.size() > 3);
            for (int i = 0; i < (int)words.size(); i++) {
                DOUBLE value = LogEnergy(atoi(words[i].c_str()));
                if (pwords[i+1] == "cu") logML_BASE = value;
                else if (pwords[i+1] == "cc") logMLclosing = value;
                else if (pwords[i+1] == "ci") logMLintern = value;
                else if (pwords[i+1] == "TerminalAU") SetAU(value);
                // else cout << "Unidentified parameter: " << words[i] << " " << pwords[i+1] << endl;
            }
        } else {
            assert((int)words.size() > 5);
            logML_BASE = LogEnergy(atoi(words[0].c_str()));
            logMLclosing = LogEnergy(atoi(words[2].c_str()));
            logMLintern = LogEnergy(atoi(words[4].c_str()));
        }
        break;
    }
}
void Convert::SetAU(DOUBLE value)
{
    logTermAU = value;
    inittermau = true;
}
void Convert::SetAU(string value)
{
    logTermAU = LogEnergy(atoi(value.c_str()));
    inittermau = true;
}
void Convert::ReadMisc(bool lxc)
{
    string str;
    vector<string> words, pwords;
    while (getline(ifs, str)) {
        if (str == "") break;
        GetWords(str, words);
        if (str.find("*") != string::npos) continue;
        if ((int)pwords.size() == (int)words.size()+2) {
            assert((int)words.size() > 2);
            for (int i = 0; i < (int)words.size(); i++) {
                if (lxc) {
                    if (pwords[i+1] == "LXC" || pwords[i+1] == "lxc") lxc37 = atof(words[i].c_str());
                } else if (pwords[i+1] == "TerminalAU") SetAU(words[i]);
            }
        } else {
            assert((int)words.size() >= 4);
            if (lxc) {
                if ((int)words.size() > 4)
                    lxc37 = atof(words[4].c_str());
            } else SetAU(words[2]);
        }
        break;
    }
}

void Convert::ReadString(DOUBLE* array, string & loopstr)
{
    string str;
    vector<string> words;
    for (int i = 0; getline(ifs, str); i++) {
        if (str == "") break;
        else if (str.find("*") != string::npos) {
            i--; continue;
        }
        GetWords(str, words);
        assert((int)words.size() > 1);
        loopstr += words[0] + " ";
        array[i] = LogEnergy(atoi(words[1].c_str()));
    }
}

void Convert::FillINF()
{
    fill(&(logstack[0][0]), &(logstack[0][0])+(NBPAIRS+1)*(NBPAIRS+1), -INF);
    fill(&(logmismatchH[0][0][0]), &(logmismatchH[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logmismatchI[0][0][0]), &(logmismatchI[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logmismatch1nI[0][0][0]), &(logmismatch1nI[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logmismatch23I[0][0][0]), &(logmismatch23I[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logmismatchM[0][0][0]), &(logmismatchM[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logmismatchExt[0][0][0]), &(logmismatchExt[0][0][0])+(NBPAIRS+1)*5*5, -INF);
    fill(&(logdangle5[0][0]), &(logdangle5[0][0])+(NBPAIRS+1)*5, -INF);
    fill(&(logdangle3[0][0]), &(logdangle3[0][0])+(NBPAIRS+1)*5, -INF);
    fill(&(logint11[0][0][0][0]), &(logint11[0][0][0][0])+(NBPAIRS+1)*(NBPAIRS+1)*5*5, -INF);
    fill(&(logint21[0][0][0][0][0]), &(logint21[0][0][0][0][0])+(NBPAIRS+1)*(NBPAIRS+1)*5*5*5, -INF);
    fill(&(logint22[0][0][0][0][0][0]), &(logint22[0][0][0][0][0][0])+(NBPAIRS+1)*(NBPAIRS+1)*5*5*5, -INF);
    fill(&(loghairpin[0]), &(loghairpin[0])+MAXLOOP+1, -INF);
    fill(&(logbulge[0]), &(logbulge[0])+MAXLOOP+1, -INF);
    fill(&(loginternal[0]), &(loginternal[0])+MAXLOOP+1, -INF);
    fill(&(logninio)[0], &(logninio[0])+MAXLOOP+1, -INF);
    fill(&(Triloop[0]), &(Triloop[0])+40, -INF);
    fill(&(Tetraloop[0]), &(Tetraloop[0])+40, -INF);
    fill(&(Hexaloop[0]), &(Hexaloop[0])+40, -INF);
}

void Convert::ReadOnlyMisc(string& file)
{
    string str;
    ifs.open(file.c_str());
    while (getline(ifs, str)) {
        int num = param_type(str);
        if (num == Misc) {
            ReadMisc(true);
            break;
        }
    }
    ifs.close();
    FillINF();
}

bool Convert::ConvertParamFile(string& file)
{
    string str;
    ReadOnlyMisc(file);
    ifs.open(file.c_str());
    old_param = true;
    if (getline(ifs, str) &&
        (str.find("## RNAfold parameter file v2.0") != string::npos || file.find("parameters_BLstar_Vienna.par") != string::npos))
        old_param = false;
    while (getline(ifs, str)) {
        int num = param_type(str);
        if (num == Stack) {
            if (old_param) Read2Dim(&(logstack[0][0]), NBPAIRS+1, NBPAIRS+1, 0, 0);
            else Read2Dim(&(logstack[0][0]), NBPAIRS+1, NBPAIRS+1, 1, 1);
        } else if (num == MisH) {
            if (old_param) Read3Dim(&(logmismatchH[0][0][0]), NBPAIRS+1, 5, 5, 0, 0, 0);
            else Read3Dim(&(logmismatchH[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == MisI) {
            if (old_param) Read3Dim(&(logmismatchI[0][0][0]), NBPAIRS+1, 5, 5, 0, 0, 0);
            else Read3Dim(&(logmismatchI[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == Mis1n) {
            Read3Dim(&(logmismatch1nI[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == MisI23) {
            Read3Dim(&(logmismatch23I[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == MisM) {
            Read3DimSmooth(&(logmismatchM[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == MisE) {
            Read3DimSmooth(&(logmismatchExt[0][0][0]), NBPAIRS+1, 5, 5, 1, 0, 0);
        } else if (num == Dan5) {
            if (old_param) Read2DimSmooth(&(logdangle5[0][0]), NBPAIRS+1, 5, 0, 0);
            else Read2DimSmooth(&(logdangle5[0][0]), NBPAIRS+1, 5, 1, 0);
        } else if (num == Dan3) {
            if (old_param) Read2DimSmooth(&(logdangle3[0][0]), NBPAIRS+1, 5, 0, 0);
            else Read2DimSmooth(&(logdangle3[0][0]), NBPAIRS+1, 5, 1, 0);
        } else if (num == Int11) {
            Read4Dim(&(logint11[0][0][0][0]), NBPAIRS+1, NBPAIRS+1, 5, 5,
                                    1, 1, 0, 0);
        } else if (num == Int21) {
            Read5Dim(&(logint21[0][0][0][0][0]), NBPAIRS+1, NBPAIRS+1, 5, 5, 5,
                                    1, 1, 0, 0, 0);
        } else if (num == Int22) {
            if (old_param) {
                Read6Dim(&(logint22[0][0][0][0][0][0]), NBPAIRS+1, NBPAIRS+1, 5, 5, 5, 5,
                                        1, 1, 1, 1, 1, 1,
                                        0, 0, 0, 0, 0, 0);
            } else {
                Read6Dim(&(logint22[0][0][0][0][0][0]), NBPAIRS+1, NBPAIRS+1, 5, 5, 5, 5,
                                        1, 1, 1, 1, 1, 1,
                                        1, 1, 0, 0, 0, 0);
            }
        } else if (num == Hairpin) {
            Read1Dim(&(loghairpin[0]), MAXLOOP+1, 0);
        } else if (num == Bulge) {
            Read1Dim(&(logbulge[0]), MAXLOOP+1, 0);
        } else if (num == Interior) {
            Read1Dim(&(loginternal[0]), MAXLOOP+1, 0);
        } else if (num == Ninio) {
            Readninio();
        } else if (num == ML) {
            ReadML();
        } else if (num == Misc) {
            ReadMisc();
        } else if (num == Tri) {
            ReadString(&(Triloop[0]), Triloops);
        } else if (num == Tetra) {
            ReadString(&(Tetraloop[0]), Tetraloops);
        } else if (num == Hexa) {
            ReadString(&(Hexaloop[0]), Hexaloops);
        } else
            ;
    }
    ifs.close();
    return true;
}

}
}
