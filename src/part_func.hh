#ifndef _PART_FUNC_HH
#define _PART_FUNC_HH

#include "param.hh"
#include "matrix.hh"
#include "plot_struct.hh"
#include "centroid.hh"
#include "douter_concat.hh"
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <iterator>
#include <iomanip>


#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


#define PREC (15)
#define JUMP (1)
#define TMP (HOMEDIR + string("outer/"))
#define STMP (HOMEDIR + string("prob/"))
#define TYPE (6)

#define BPPM(pos, dist, start) (bppm[(pos)-(start)][(dist)])
#define PROM(pos, type, start) (prom[((pos)-(start))*TYPE+type])

namespace Rfold {

using namespace Parameter;
using std::ostringstream;
using std::istringstream;
using std::istream_iterator;
using std::istreambuf_iterator;
using std::ofstream;
using std::reverse;
using std::min;
using std::max;
using std::abs;
using std::setprecision;
using std::setiosflags;
using std::transform;
using std::back_inserter;
using std::max_element;
using std::numeric_limits;
using std::accumulate;
using std::find_if;
using std::distance;
using std::plus;



class Arg {
public:
    Arg() {
        str = "";
        input = "";
        outer_input = "";
        output = "";
        name = "";
        ene = "";
        mout = "";
        hard_const = "";
        start = -1;
        end = -1;
        length = 0;
        seqID = 0;
        id = -1;
        chunk = 2;
        constraint = 200;
        window = 10;
        mpoint = -1;
        mtype = -1;
        init_calc = Calc::Pre;
        temp = temperature;
        gamma = 1.0;
        minp = -1.0;
        eSDC_flag = false;
        stem_flag = false;
        acc_flag = false;
        prof_flag = false;
        debug = false;
        keep_flag = false;
        stemdb_flag = false;
        compl_flag = false;
        mea_flag = false;
        pre_flag = false;
        save_memory = false;
        init_file = false;
        text = false;
        outtext = false;
        image = false;
        mout_flag = false;
        cd = false;
        boundary = false;
        entro_flag = false;
        mfe_flag = false;
    }
    virtual ~Arg() {}
    enum Calc {Divide, Connect, Stemdb, Bpp, Pre,};
    enum Mut {Sub, Del, Ins};
    string str;
    string input;
    string outer_input;
    string output;
    string name;
    string ene;
    string mout;
    string hard_const;
    LEN start;
    LEN end;
    LEN length;
    int seqID;
    int id;
    int chunk;
    int constraint;
    int window;
    int mpoint;
    int mtype;
    int init_calc;
    int temp;
    DOUBLE gamma;
    DOUBLE minp; // minimum base pairing probability to output
    bool eSDC_flag;
    bool stem_flag;
    bool acc_flag;
    bool prof_flag;
    bool keep_flag; // when connecting, keep temp file;
    bool stemdb_flag; // when connecting, output stem probability db;
    bool compl_flag;  // calculate complementary sequence;
    bool mea_flag;
    bool pre_flag;
    bool save_memory;
    bool init_file;
    bool debug;
    bool text;
    bool outtext;
    bool image;
    bool mout_flag;
    bool cd;
    bool boundary;
    bool entro_flag;
    bool mfe_flag;

};

class ParasoR {
private:
    enum Out {STEM, ACC, PROF, BPPIM, BPP, MOTIF, PROFIM, ENTRO};
    enum File {Part, Shrunk, Whole};
    DOUBLE gamma;

    /* part_func.cc */

    DOUBLE GetInStem(LEN, LEN);
    DOUBLE GetInMultiBif(LEN, LEN);
    DOUBLE GetInMulti2(LEN, LEN);
    DOUBLE GetInMulti1(LEN, LEN);
    DOUBLE GetInMulti(LEN, LEN);
    DOUBLE GetInStemend(LEN, LEN);
    void SetInsideMat(LEN, LEN);
    // void CalcInside();

    DOUBLE GetOutStem(LEN, LEN);
    DOUBLE GetOutMultiBif(LEN, LEN);
    DOUBLE GetOutMulti2(LEN, LEN);
    DOUBLE GetOutMulti1(LEN, LEN);
    DOUBLE GetOutMulti(LEN, LEN);
    DOUBLE GetOutStemend(LEN, LEN);
    // void SetOutsideMat(LEN, LEN);

    DOUBLE GetStemDelta(LEN, LEN, bool);
    DOUBLE GetStemDelta(LEN, LEN, int, bool);
    DOUBLE CalcDalpha(LEN, LEN, int, DOUBLE);
    DOUBLE CalcDbeta(LEN, LEN, int, DOUBLE);
    void CalcInDeltaOuter(LEN);
    void CalcOutDeltaOuter(LEN);
    void CalcChunkInside(bool = true);
    void CalcChunkOutside();
    void CalcDeltaInOut(bool);
    void CalcDeltaInOut();

    void ReadVec(Vec&, string&);
    bool ReadBinVec(int, Vec&, ifstream&);
    bool ReadDouterInside(Mat&, string, Vec&);
    bool ReadBinDouterInside(Mat&, string, Vec&);
    bool ReadDouterOutside(Mat&, string, Vec&);
    bool ReadBinDouterOutside(Mat&, string, Vec&);
    void GetSumDouterList(const Vec& old_douter, Vec& sum_douter, bool);
    DOUBLE GetDenominator(const Vec& douter, const Vec& sum_douter, bool);
    DOUBLE GetNumerator(const Mat& douter, const Vec& sum_douter, LEN k, bool inside);
    void ConnectInDo(Vec&, Mat&, int, string, bool, bool = true);
    void ConnectOutDo(Vec&, Mat&, int, string, bool, bool = true);
    bool ConnectDo(bool, bool);
    void ConnectDo(bool);
    void ConnectDoSaved(bool);
    void ReadBinFirstDouter(bool, Vec&, string&);
    void ReadFirstDouter(bool, Vec&, string&);
    void ReadFirstDouter(bool);
    void ConnectDoSaved(bool, Vec&);
    bool ConnectSavedFiles(bool, bool);
    bool ConnectSavedFiles(bool = false);

    void WriteBpp(Mat&);
    void WriteAcc(Mat&);
    void WriteStemProb(Vec&, const Mat&);
    void WriteStemProb(Vec&);
    void OutputStemProb(bool);
    void OutputBppCond(bool);

    string GetDoFile(bool);
    string GetStemFile(bool, bool = false, bool = false);
    string GetDividedStemFile(bool, bool = false, bool = false);
    string GetIDFileList(string, int tid);
    string GetShrunkFileList(int, bool, int tid = -1);
    string GetTempFileList(bool, int tid = -1);
    void WriteDouterTemp(ofstream&, Mat&);
    void WriteBinDouterTemp(ofstream&, Mat&);
    void StoreDouterTemp(bool);
    void LeaveDouterEdgeRegion(bool);
    void StoreDouterTempShrunk(bool);
    void StoreDouter(string, Vec&, bool, bool);
    void StoreStem(string, Vec&, bool);
    void StoreProf(string, vector<char>&);

    void PrintMat(bool);
    void InitMatrix(LEN, LEN);
    void Init(bool full = false, bool connect = false);
    void InitBpp(bool full = false, bool set = true);
    void SetSequence(const string&, int, LEN);
    void SetSequence(const string&, bool = true);
    void RemoveTempFiles(bool, int = 0);
    void RemoveStem(bool, bool = false);
    void ConcatDo();
    static void Connect(ParasoR&, bool, bool = false);

    /* part_func_prof */

    DOUBLE HairpinAcc(LEN, LEN);
    DOUBLE InnerLoopAcc(LEN, LEN, LEN, LEN, bool = false);
    DOUBLE InteriorAcc(LEN, LEN, bool = true, bool = false);
    DOUBLE BulgeAcc(LEN, LEN);
    DOUBLE MultiAcc(LEN, LEN);
    DOUBLE Multi2Acc(LEN, LEN);
    DOUBLE Otherprof(Vec&);
    DOUBLE Stemprof(LEN);
    DOUBLE Bulgeprof(LEN);
    DOUBLE Interiorprof(LEN);
    DOUBLE Hairpinprof(LEN);
    DOUBLE Multiprof(LEN);
    DOUBLE InteriorLinearLeft(LEN x1, LEN x2, LEN i, LEN k, bool exbulge, bool exinter);
    DOUBLE InteriorLinearRight(LEN x1, LEN x2, LEN j, LEN l, bool exbulge, bool exinter);
    DOUBLE InteriorLinear(LEN, LEN, bool, bool);
    DOUBLE StemLinear(LEN, LEN);
    DOUBLE BulgeLinear(LEN, LEN);
    DOUBLE InteriorLinear(LEN, LEN);
    DOUBLE HairpinLinear(LEN, LEN);
    DOUBLE MultiLinear(LEN);
    void AddCumProf(Vec& cump, LEN pos, Vec& prof);
    void GetProfsLinearRight(LEN pos, Vec& prof, LEN bstart, LEN end);
    void GetProfsLinear(LEN, Vec&, LEN, LEN, DOUBLE);
    void GetProfs(LEN, Vec&, LEN = 1, DOUBLE = 0., bool = false);

    /* part_func_delta */

    DOUBLE GetPairingStruct(LEN, LEN, DOUBLE, DOUBLE);
    DOUBLE ReCalcDouterInside(LEN);
    DOUBLE ReCalcDouterOutside(LEN);

    void PreCalcInside(LEN);
    void PreCalcOutside(LEN);
    void PreCalcOutside(LEN, DOUBLE);

    DOUBLE ReproduceLocalOuter(LEN);
    DOUBLE ExpandLocalOuter(DOUBLE, LEN, LEN, LEN);
    DOUBLE ShrinkLocalOuterInside(DOUBLE, LEN);
    DOUBLE ShrinkLocalOuter(DOUBLE, LEN);
    DOUBLE RewindLocalOuter(LEN, DOUBLE);
    DOUBLE SlideLocalOuter(LEN, DOUBLE);

    DOUBLE GetOutStem(LEN, LEN, DOUBLE);
    void SetOutsideMat(LEN, LEN, DOUBLE);
    int ReadPartConnectedDouter(bool, Vec&);
    int ReadBinPartConnectedDouter(bool, Vec&);
    bool ReadConnectedDouter(bool);
    bool ReadBinConnectedDouter(bool);
    void ReadStem(Vec&, string, LEN, LEN);
    void ReadBinStem(Vec&, string, LEN, LEN);
    void ReadStemVec(Vec&, bool);

    void InitRowMat(Matrix&, LEN);
    void InitColMat(Matrix&, LEN);

    void CalcInside(LEN);
    void CalcInsideFromRight(LEN);
    void CalcOutside(LEN);
    void CalcOutside(LEN, DOUBLE);
    void CalcForward(LEN);
    void CalcForward(LEN, DOUBLE);
    void SetRawRangedMatrix(bool);
    DOUBLE SetRangedMatrix(LEN, bool set = false);
    void StoreBppSlide(LEN, LEN, Mat&);
    void StoreBppSlide(LEN, LEN, Vec&);
    void StoreAccProfSlide(LEN, LEN, int, bool = false);
    void StoreAreaBppSlide(LEN, LEN, Mat&);
    void StoreAreaBppSlide(LEN, LEN, Vec&);
    void CheckDouter(LEN, LEN, DOUBLE);
    void SetProbs(Mat&);
    void SetProbs(Vec&);
    void CalcSlidingWindowBoundary(Vec&, LEN, LEN, bool set = true);
    template <class Probs>
    void CalcSlidingWindowStem(Probs&, LEN, LEN, bool set = true);
    void CalcSlidingWindowAcc(Vec&, int, LEN, LEN, bool set = true);

    template <class Substructure>
    void CalcSlidingWindowProf(Substructure&, LEN, LEN, bool set = true);
    void CalcStem(Mat&);
    void CalcStem(Vec&, bool = false);
    void CalcStem(bool = false, bool = false);
    void CalcAcc(Vec&, int);
    void CalcAcc(bool = false);
    void CalcProf(Vec& P);
    void CalcProf(vector<char>& P);
    void CalcProf(bool = false, bool = false);
    void DrawImage(vector<int>&, Vec&, string&, int = 1);
    bool GetWholeImage(string, vector<int>&, Vec&, int = 1);
    void GetImage(string, LEN, LEN, Vec, int = 1);
    void CalcBpp(bool = false, DOUBLE = -1.0, bool = false);
    void CalcEntropy(bool = true);
    void CalcMEA(bool = true, bool = false, bool = false, bool = false);
    void CalcRangeStem(Vec&, LEN, LEN, int);
    void CalcRangeAcc(Vec&, int, LEN, LEN, int);

    void TabToNewLine(string&);
    bool ReadStemToSingleFile(string&, string&, bool, string);
    bool ReadBinStemToSingleFile(string&, string&, bool);
    void ConcatStemdb(bool, bool = false);

    void CalcInsideOuter(LEN);
    void CalcOutsideOuter(LEN);
    void CalcOuter();
    void CalcBppAtOnce(int out, bool image, DOUBLE thres);
    void CalcAllAtOnce(int, bool = false, DOUBLE = 0.0, bool = false);
    void PrintRfoldHeader(int);

    void SetOriginalDouter(Vec&, Vec&);
    void CopyOriginalOuter(Vec&, Vec&);
    void SlidingMatrixForMut(LEN, int);
    void ComeBackChange(LEN, int, int);
    bool SlidingSequence(LEN, int, int);
    void Recalculation(LEN);
    DOUBLE MaxDiff(const Vec&, const Vec&);
    // DOUBLE Correlation(const Vec&, const Vec&);
    string GetDiffFile(int, bool);
    void CutOffVector(int, Vec&, LEN);
    void EditVector(int, Vec&, Vec&, LEN);
    Vec GetSubstrVec(const Vec&, LEN, int);
    void WriteDiffFile(LEN, int, DOUBLE, DOUBLE, DOUBLE, string&, bool, bool);
    void WriteDiff(int, LEN, int, Vec, Vec&, bool, bool, bool, string&, int);
    void ChangeBase(LEN, Arg&, bool&);
    void MutatedStem(Arg&);

    /* part_func_hc.cc */
    DOUBLE GetInStemHC(LEN, LEN);
    DOUBLE GetInMultiBifHC(LEN, LEN);
    DOUBLE GetInMulti2HC(LEN, LEN);
    DOUBLE GetInMulti1HC(LEN, LEN);
    DOUBLE GetInMultiHC(LEN, LEN);
    DOUBLE GetInStemendHC(LEN, LEN);

    DOUBLE GetOutStemHC(LEN, LEN);
    DOUBLE GetOutStemHC(LEN, LEN, DOUBLE);
    DOUBLE GetOutMultiBifHC(LEN, LEN);
    DOUBLE GetOutMulti2HC(LEN, LEN);
    DOUBLE GetOutMulti1HC(LEN, LEN);
    DOUBLE GetOutMultiHC(LEN, LEN);
    DOUBLE GetOutStemendHC(LEN, LEN);

    void CalcInsideOuterHC(LEN);
    void CalcOutsideOuterHC(LEN);
    void CalcOuterHC();
    bool SameSequenceHC(LEN, LEN);
    bool InnerLoopHC(LEN, LEN, LEN, LEN);
    bool IsPairHC(LEN, LEN);
    bool IsLoopHC(LEN);

    LEN RightRange(LEN i) {
        return min(seq.length, i+_constraint+1);
    }
    LEN RightBpRange(LEN i) {
        return min(seq.length, i+_constraint);
    }
    LEN LeftRange(LEN j) {
        return max(static_cast<LEN>(0), j-_constraint-1);
    }
    LEN LeftBpRange(LEN j) {
        // return max(static_cast<LEN>(0, j-_constraint);
        return max(static_cast<LEN>(1), j-_constraint);
    }
    bool IsRange(LEN i, LEN j) {
        return (i >= j-_constraint-1 && i > 0 && j < seq.length);
    }
    bool IsOnlyRange(LEN i, LEN j) {
        return (i >= j-_constraint-1 && i >= 0 && j <= seq.length);
    }

    int BestChunk(int tchunk) {
        if (seq.length/static_cast<LEN>(tchunk) <= 3*_constraint) return seq.length/(_constraint*3)+1;
        else return tchunk;
    }
    int GetColumn(ifstream& ifs) {
        int temp = 0;
        ifs.read((char*)&temp, sizeof(int));
        return temp;
    }
public:
    ParasoR() {
        id = 0;
        chunk = 0;
        InitEnergyParam();
        cut = false;
        delta = true;
        binary = true;
        memory = false;
        mfe = false;
        hard = false;
        min_or_sum = Logsumexp;
    }
    ParasoR(string& name) : name(name) {
        id = 0;
        chunk = 0;
        InitEnergyParam();
        cut = false;
        delta = true;
        binary = true;
        memory = false;
        mfe = false;
        hard = false;
        min_or_sum = Logsumexp;
    }
    virtual ~ParasoR() {}
    int id;
    int chunk;
    LEN _length;    // the length of original sequence;
    LEN _constraint;
    int _window;    // window := average width, region := window-1;
    LEN _start;
    LEN _end;
    string name;
    Sequence seq;
    Matrix alpha;
    Matrix beta;
    Vec bppv;       // used for original stem probabilities or accessibilities;
    Mat bppm;       // used in stem probability calculation;
    Vec prom;       // used in profile calculation;
    bool cut;       // cut a needless region of sequence;
    bool delta;     // douter or outer flag;
    bool binary;    // binary storage flag;
    bool memory;    // memory saving mode;
    bool mfe;       // output MFE structure;
    bool hard;      // hard constraint flag;
    DOUBLE (*min_or_sum)(DOUBLE, DOUBLE);
    static const bool ene = true;    // output an energy of accessibility;
    static const bool linear = true; // linear profile calculation;
    static const bool noout = false;
    static const bool print = (false && !noout);
    static const bool debug = (false && print);
    static const bool ddebug = false;
    static const LEN CONST = 4;

    DOUBLE bpp_func(LEN i, LEN j, bool deb = false) {
        if (hard) return bppHard(i, j, deb);
        else if (delta) return bppDelta(i, j, deb);
        else    return bpp(i, j, deb);
        return bpp(i, j, deb);
    }
    DOUBLE bpp(LEN, LEN, bool deb = false);
    DOUBLE bppDelta(LEN, LEN, bool deb = false);
    DOUBLE bppHard(LEN, LEN, bool deb = false);

    DOUBLE acc(LEN, LEN);
    DOUBLE accDelta(LEN, LEN, DOUBLE);
    void profile(LEN, Vec&);
    void profileDelta(LEN, DOUBLE, Vec&);
    void profileDeltaLinear(LEN, LEN, DOUBLE, Vec&);
    void RemoveTempFilesSaveMem(bool);
    static void DivideChunk(Arg&, bool shrink = true);
    static void Connect(Arg&, bool shrink = true);
    static void Stemdb(Arg&, bool shrink = true);
    static void PreviousCalculation(Arg&, bool shrink = true);
    static void main(Arg&, bool shrink = true);
    void SetText(bool text) {
        if (text) binary = false;
    }
    void SetNoDelta() {
        delta = false;
    }
    void SetMFEFlag(bool tmfe) {
        mfe = tmfe;
        if (tmfe) {
            min_or_sum = Logsummax;
        } else {
            min_or_sum = Logsumexp;
        }
    }
    void SetHardConstFlag(bool thc) {
        if (thc && seq.hclength == seq.length) {
            hard = true;
        } else {
            hard = false;
        }
    }
    void SetGamma(DOUBLE tgamma) {
        gamma = tgamma;
    }
    void SetMemory(bool tmemory) {
        memory = tmemory;
        if (static_cast<LEN>(_constraint) >= seq.length-1 || 3*static_cast<LEN>(_constraint) < seq.length) {
            memory = false;
        }
        cut = (cut || tmemory);
    }
    void SetWindow(int window, bool acc = false)
    {
        if (!noout && acc)
            cout << "# Window : " << window << endl;
        _window = window;
    }
    void SetConstraint(int constraint, LEN length)
    {
        _length = length;
        _constraint = (constraint > 0) ? min(static_cast<LEN>(constraint), _length-1) : _length-1;
    }
    void SetRawCons(int constraint, LEN length)
    {
        _length = length;
        _constraint = static_cast<LEN>(constraint);
    }
    void SetBasicParam(int constraint, LEN length, string& tname, bool shrink)
    {
        name = tname;
        if (shrink) SetConstraint(constraint, length);
        else SetRawCons(constraint, length);
    }
    void SetBasicParam(Arg& arg, bool shrink = false)
    {
        SetBasicParam(arg.constraint, max(static_cast<LEN>(arg.str.length()), arg.length), arg.name, shrink);
        if (static_cast<LEN>(arg.str.length()) > 0) {
            SetSequence(arg.str);
            if (static_cast<LEN>(arg.hard_const.length()) == seq.length) {
                seq.SetHardConstraint(arg.hard_const);
            } else if (arg.hard_const.length() > 0) {
                if (!noout) cout << "# Hard const length error. " << static_cast<LEN>(arg.hard_const.length()) << " != " << seq.length << endl;
            }
        } else {
            SetSequence(arg.input, arg.seqID, arg.length);
        }
        SetText(arg.text);
        SetGamma(arg.gamma);
        SetMemory(arg.save_memory);
    }
    void SetRange(LEN start, LEN end) {
        _start = max(static_cast<LEN>(1), start);
        _end = min(seq.length, static_cast<LEN>(end));
    }
    void SetBpRange(LEN start, LEN end) {
        _start = max(static_cast<LEN>(0), start);
        _end = min(seq.length, static_cast<LEN>(end));
    }
    bool SetChunkId(int tid, int tchunk, bool range = true)
    {
        tchunk = min(tchunk, BestChunk(tchunk));
        if (tid >= tchunk) return false;
        else {
            id = tid;
            chunk = tchunk;
            if (range) {
                _start = (seq.length/tchunk) * (tid);
                _end = (tid+1 == tchunk) ? seq.length : (seq.length/tchunk) * (tid+1);
                if (!noout) cout << "#--(s, t) " << _start << " " << _end << endl;
            }
            Init();
            return true;
        }
    }
    void SetIndex(LEN start, LEN end, bool delta, bool set = true)
    {
        alpha.SetIndex(start, end, delta, set);
        beta.SetIndex(start, end, delta, set);
    }
    void InitVec(Mat& data, LEN pos) {
        data[pos].assign(data[pos].size(), -INF);
    }
    static DOUBLE expInf(DOUBLE value) {
        if (Is_INF(value))  return 0.;
        return exp(value);
    }
    static DOUBLE ProbEntropy(DOUBLE p) {
        if (p == static_cast<DOUBLE>(0.0)) return static_cast<DOUBLE>(0.);
        return p*log2(p);
    }
    static char base(int c) {
        switch(c) {
            case 1: return 'A';
            case 2: return 'C';
            case 3: return 'G';
            case 4: return 'U';
            default: return 'N';
        }
    }
    static char substructure(int c) {
        switch(c) {
            case 0: return 'B';
            case 1: return 'O';
            case 2: return 'H';
            case 3: return 'M';
            case 4: return 'S';
            default: return 'I';
        }
    }
    static void AppendProf(vector<char>& P, Vec& prof) {
        for (Vec::iterator it = prof.begin(); it != prof.end(); it += TYPE) {
            P.push_back(substructure(std::distance(it, max_element(it, it+TYPE))));
        }
    }
    static void AppendProf(Vec& P, Vec& prof) {
        P = prof;
        // for (Vec::iterator it = prof.begin(); it != prof.end(); it++)
        //     P.push_back(*it);
    }

};

}

#endif
