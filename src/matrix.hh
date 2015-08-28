
#ifndef _MATRIX_HH
#define _MATRIX_HH

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <iterator>
#include "energy_const.hh"

#define Ind(mat, j) ((j)%(mat.length+1))

#define Stem(mat, i, j) (mat.stem[Ind(mat, j)][(j)-(i)])
#define Stemend(mat, i, j) (mat.stemend[Ind(mat, j)][(j)-(i)])
#define Multi(mat, i, j) (mat.multi[Ind(mat, j)][(j)-(i)])
#define Multi1(mat, i, j) (mat.multi1[Ind(mat, j)][(j)-(i)])
#define Multi2(mat, i, j) (mat.multi2[Ind(mat, j)][(j)-(i)])
#define Multibif(mat, i, j) (mat.multibif[Ind(mat, j)][(j)-(i)])
#define GetDo(mat, j, h) ((j-mat.istart < 0) ? (0.0) : ((j > mat.iend) ? (0.0) : (mat.douter[j-mat.istart][h])))
#define DOuter(mat, j, h) (mat.douter[j-mat.istart][h])
#define DOuterInd(mat, j) (j-mat.istart)
#define Outer(mat, j) (mat.outer[j-mat.istart])



namespace Rfold {

using std::string;
using std::vector;
using std::deque;
using std::min;
using std::max;
using std::transform;
using std::ifstream;
using std::istringstream;
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::back_inserter;
using std::distance;
using std::swap;


// typedef vector<DOUBLE> Vec;
// typedef vector<Vec> Mat;

using Vec = vector<DOUBLE>;
using Mat = vector<Vec>;

struct Base {
    char operator()(char c) {
        switch(c) {
            case 'A': case 'a': return 1;
            case 'C': case 'c': return 2;
            case 'G': case 'g': return 3;
            case 'T': case 'U': case 't': case 'u': return 4;
            case 'X': return 5;
            case 'K': return 6;
            case 'I': return 7;
            case '$': case 'N': return 0;
            default: return 0;
        }
    }
};


template<class Data>
void PrintVec(const vector<Data>& elem, bool end = true)
{
    ostream_iterator<Data> out_it(cout, ", ");
    if ((LEN)elem.size() > 1)
        copy(elem.begin(), elem.begin()+(LEN)elem.size()-1, out_it);
    if ((LEN)elem.size() > 0)
        cout << *(elem.begin()+(LEN)elem.size()-1);
    if (end) cout << endl;
}

class Matrix {
private:
    void Initialize();
public:
    static const bool debug = false;
    bool _innerset;
    bool _inside;
    int _constraint; // height;
    LEN length;     // width;
    LEN olength;    // width of outer;
    LEN istart;      // index of start position on full sequence;
    LEN iend;        // index of end position on full sequence;

    enum Type { Stem, Stemend, Multi, Multi1, Multi2, Multibif };
    Vec outer;
    Mat douter;
    Mat stem;
    Mat stemend;
    Mat multi;
    Mat multi1;
    Mat multi2;
    Mat multibif;

    Matrix() {
       _innerset = false;
    }
    Matrix(LEN length, int constraint, bool inside)
     : _inside(inside), _constraint(constraint), length(length) {
        Initialize();
    }
    Matrix(LEN length, int constraint, bool inside, Vec& touter) {
        Initialize();
        outer = touter;
    }
    virtual ~Matrix(){}
    void operator=(const Matrix& right) {
        outer = right.outer;
        douter = right.douter;
        stem = right.stem;
        stemend = right.stemend;
        multi = right.multi;
        multi1 = right.multi1;
        multi2 = right.multi2;
        multibif = right.multibif;
        length = right.length;
        _constraint = right._constraint;
        _inside = right._inside;
        _innerset = right._innerset;
    }
    static void PrintMat(const Mat&, const string&);
    void Insert(LEN x) {
        outer.insert(outer.begin()+x-istart, 0.0);
        iend++;
        olength++;
    }
    void Delete(LEN x) {
        outer.erase(outer.begin()+x-istart);
        iend--;
        olength--;
    }
    void Print(const string&);
    void SetIndex(LEN s, LEN e, bool delta = true, bool set = true) {
        istart = s;
        iend = e;
        olength = iend-istart;
        if (_inside && debug)
            cerr << "start, end: " << s << " " << e << endl;
        if (delta) {
            douter = Mat(olength+1, Vec(_constraint+1, 0));
        }  else if (set) {
            outer = Vec(olength+1, 0);
        }
    }
    bool isSet() {
        return _innerset;
    }

};
}


#endif