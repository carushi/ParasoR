#ifndef _CENTROID_HH
#define _CENTROID_HH

#include "matrix.hh"

namespace Centroid {

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::max;
using std::abs;

/**
 * Gets structure string in bracket style.
 * @param bppm base pairing probability matrix.
 * @param end the last index of of included bases.
 * @param gamma centroid parameter.
 * @return structure string
 */
static string GetMEAStructure(Rfold::Mat& bppm, LEN end, DOUBLE gamma)
{
    string str;
    str.assign(end, '.');
    for (LEN j = 0; j < bppm.size(); j++) {
        for (LEN dist = 0; dist < bppm[j].size(); dist++) {
            if (bppm[j][dist] > 1./(gamma+1.)) {
                if (j-dist >= 0)
                    str[j-dist] = '(';
                if (j < end)
                    str[j] = ')';
            }
        }
    }
    return str;
}

static bool IsSame(DOUBLE a, DOUBLE b) {
    const DOUBLE EPS = 0.000001;
    return abs(a-b) < EPS;
}

/**
 * Set () to str at centroid base pair position.
 * @param i, j range.
 * @param start, end sequence start and end.
 * @param m dp table.
 * @param bppm bpp table.
 * @param gamma centroid parameter.
 */
static void SetStructure(LEN i, LEN j, LEN start, LEN end, Rfold::Mat& m, string& str, Rfold::Mat& bppm, DOUBLE gamma)
{
    while (TURN < j-i) {
        if (IsSame(m[i][j], m[i+1][j])) i++;
        else if (IsSame(m[i][j], m[i][j-1])) j--;
        else if (i+1 < j-1 && IsSame(m[i][j], m[i+1][j-1]+(gamma+1.)*bppm[j-1][j-i]-1.)) {
            str[i-1] = '('; str[j-1] = ')';
            i++; j--;
        } else {
            for (LEN k = i+1; k+1 < j; k++) {
                if (IsSame(m[i][j], m[i][k]+m[k+1][j])) {
                    SetStructure(i, k, start, end, m, str, bppm, gamma);
                    SetStructure(k+1, j, start, end, m, str, bppm, gamma);
                    return;
                }
            }
            cout << "Error: No match during traceback." << endl;
            exit(1);
        }
    }
}

/**
 * Centroid structure prediction.
 * @param start, end sequence start and end.
 * @param bppm bpp table.
 * @param gamma centroid parameter.
 * @return gamma centroid structure.
 */
static string GetCentroidStructure(Rfold::Mat& bppm, LEN start, LEN end, DOUBLE gamma)
{
    start = max(static_cast<LEN>(1), start);
    LEN constraint = bppm[0].size();
    string str;
    str.assign(end, '.');
    Rfold::Mat m = Rfold::Mat(end+1, Rfold::Vec(end+1, 0));
    // for (LEN i = 0; i < bppm.size(); i++)
    //     Rfold::PrintVec(bppm[i]);
    for (LEN j = start; j <= end; j++) {
        for (LEN i = j-1; i >= 1; i--) {
            if (i+1 < j) m[i][j] = max(m[i][j], m[i+1][j]);
            if (j-1 > i) m[i][j] = max(m[i][j], m[i][j-1]);
            if (j-i > TURN && j-i <= constraint && i > 0 && j <= end) {
                if (i < j) {
                    m[i][j] = max(m[i][j], m[i+1][j-1]+(gamma+static_cast<DOUBLE>(1.))*bppm[j-1][j-i]-static_cast<DOUBLE>(1.));
                } else {
                    m[i][j] = max(m[i][j], m[i+1][j-1]+(gamma+static_cast<DOUBLE>(1.))*bppm[i-1][i-j]-static_cast<DOUBLE>(1.));
                }
            }
            for (LEN k = i+1; k < j; k++) {
                m[i][j] = max(m[i][j], m[i][k]+m[k+1][j]);
            }
        }
    }
    // for (LEN i = 0; i < m.size(); i++)
    //     Rfold::PrintVec(m[i]);
    SetStructure(start, end, start, end, m, str, bppm, gamma);
    return str;
}

/**
 * Produces stem probability vector from base pairing probability matrix.
 * @param bppm base pairing probability matrix.
 * @param end the last index of of included bases.
 * @param stem empty vector for storing stem probability.
 */
static void GetStemProb(Rfold::Mat& bppm, LEN end, Rfold::Vec& stem)
{
    stem = Rfold::Vec(end, 0);
    for (LEN j = 0; j < bppm.size(); j++) {
        for (LEN dist = 0; dist < bppm[j].size(); dist++) {
            if (bppm[j][dist] < 0)    continue;
            if (j-dist >= 0) stem[j-dist] += bppm[j][dist];
            if (j < end) stem[j] += bppm[j][dist];
        }
    }
}

/**
 * Gets base pairing probability of MEA structure.
 * @param bppm base pairing probability matrix.
 * @param end the last index of of included bases.
 * @param gamma centroid parameter.
 * @param cbpp empty vector for storing base pairing probability of MEA structure.
 */
static void GetMEABpp(Rfold::Mat& bppm, LEN end, DOUBLE gamma, vector<int>& cbpp)
{
    cbpp = vector<int>(end, 0);
    for (LEN i = 0; i < end; i++) cbpp[i] = i;
    for (LEN j = 0; j < bppm.size(); j++) {
        for (LEN dist = 0; dist < bppm[j].size(); dist++) {
            if (bppm[j][dist] > 1./(gamma+1.)) {
                if (j-dist >= 0) cbpp[j-dist] = j;
                if (j < end) cbpp[j] = j-dist;
            }
        }
    }
}

/**
 * Gets base pairing probability of string.
 * @param str structure in bracket style.
 * @param cbpp empty vector for storing base pairing probability of MEA structure.
 */
static void GetStringBpp(string& str, vector<int>& cbpp)
{
    cbpp = vector<int>(str.length(), 0);
    for (int i = 0; i < str.length(); i++) cbpp[i] = i;
    vector<int> stack;
    for (int i = 0; i < str.length(); i++) {
        if (str[i] == '(') stack.push_back(i);
        else if (str[i] == ')') {
            if (stack.size() == 0) exit(0);
            cbpp[i] = stack[stack.size()-1];
            cbpp[stack[stack.size()-1]] = i;
            stack.erase(stack.begin()+stack.size()-1);
        }
    }
}

}

#endif
