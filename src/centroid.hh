#ifndef _CENTROID_HH
#define _CENTROID_HH

#include "matrix.hh"

namespace Centroid {

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::max;

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

}

#endif