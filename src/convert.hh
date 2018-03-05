#ifndef _CONVERT_HH
#define _CONVERT_HH

#include <algorithm>
#include <iterator>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <cassert>
#include "num_operator.hh"
#include "pair_mat.hh"


namespace Rfold {
namespace Parameter {

using std::string;
using std::ifstream;
using std::vector;
using std::istream_iterator;
using std::istringstream;
using std::copy;
using std::min;
using std::cout;
using std::endl;
using std::fill;

/**
* A Convert class.
From convert_epars.c, convert_epars.h in Vienna RNA package 2.0.7.
*/
class Convert {
private:
	ifstream ifs;
	void GetWords(string&, vector<string>&);
	int param_type(string&);
	void GetArray(DOUBLE *array, int size, bool);
	void Read1Dim(DOUBLE *array, int dim, int shift, int post);
	void Read2Dim(DOUBLE* array, int dim1, int dim2, int shift1, int shift2,
	                             int post1, int post2);
	void Read3Dim(DOUBLE *array, int dim1, int dim2, int dim3, int shift1, int shift2, int shift3,
                                 int post1, int post2, int post3);
	void Read4Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4,
                                 int shift1, int shift2, int shift3, int shift4,
                                 int post1, int post2, int post3, int post4);
	void Read5Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4, int dim5,
                                 int shift1, int shift2, int shift3, int shift4, int shift5,
                                 int post1, int post2, int post3, int post4, int post5);
	void Read6Dim(DOUBLE *array, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6,
                                 int shift1, int shift2, int shift3, int shift4, int shift5, int shift6,
                                 int post1, int post2, int post3, int post4, int post5, int post6);
	void Readninio();
	void ReadML();
	void ReadMisc(bool lxc = false);
	void ReadString(DOUBLE* array, string & str);
	void Read2DimSmooth(DOUBLE* array, int dim1, int dim2, int shift1, int shift2,
								 int post1, int post2);
	void Read3DimSmooth(DOUBLE *array, int dim1, int dim2, int dim3, int shift1, int shift2, int shift3,
                                 int post1, int post2, int post3);
	void FillINF();
	void ReadOnlyMisc(string&);
	void SetAU(DOUBLE);
	void SetAU(string);
public:
	enum Params {Stack, MisH, MisI, Mis1n, MisI23, MisM, MisE, Dan5, Dan3, Int11, Int21, Int22, Hairpin, Bulge, Interior, Ninio, ML, Misc, Tri, Tetra, Hexa, };
	Convert() {}
	virtual ~Convert() {}
	bool ConvertParamFile(string&);
};

}
}


#endif
