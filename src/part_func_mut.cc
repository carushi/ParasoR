#include "part_func.hh"

namespace Rfold {

void ParasoR::SetOriginalDouter(Vec& adouter, Vec& bdouter)
{
    if (adouter.size() == 0 || bdouter.size() == 0) {
        SetRawRangedMatrix(true);
        adouter = alpha.outer;
        bdouter = beta.outer;
    }
}

void ParasoR::CopyOriginalOuter(Vec& adouter, Vec& bdouter)
{
    alpha.outer = adouter;
    alpha.iend = alpha.istart+alpha.outer.size()-1;
    alpha.olength = alpha.iend-alpha.istart;
    beta.outer = bdouter;
    beta.iend = beta.istart+beta.outer.size()-1;
    beta.olength = beta.iend-beta.istart;
}

void ParasoR::SlidingMatrixForMut(LEN x, int mtype)
{
    if (mtype <= 0)
        ;
    else if (mtype == Arg::Mut::Del) {
        alpha.Delete(x);
        beta.Delete(x);
    } else if (mtype == Arg::Mut::Ins) {
        alpha.Insert(x);
        beta.Insert(x);
    }
}

void ParasoR::ComeBackChange(LEN x, int mtype, int c)
{
    if (mtype < 0) return;
    else if (mtype == Arg::Mut::Sub) {
        seq.Sub(x, base(c), (char)c);
    } else if (mtype == Arg::Mut::Del) {
        seq.Insert(x, base(c), (char)c);
        // alpha.Insert(x);
        // beta.Insert(x);
    } else if (mtype == Arg::Mut::Ins) {
        seq.Delete(x);
        // alpha.Delete(x);
        // beta.Delete(x);
    }
}

bool ParasoR::SlidingSequence(LEN x, int c, int mtype)
{
    if (mtype < 0) return false;
    else if (mtype == Arg::Mut::Sub) {
        return seq.Sub(x, base(c), c);
    } else if (mtype == Arg::Mut::Del) {
        // alpha.Delete(x);
        // beta.Delete(x);
        return seq.Delete(x);
    } else if (mtype == Arg::Mut::Ins) {
        // alpha.Insert(x);
        // beta.Insert(x);
        return seq.Insert(x, base(c), c);
    }
    return false;
}


void ParasoR::Recalculation(LEN x)
{
    for (LEN pos = LeftRange(x); pos <= min(_end+_constraint*CONST-1, seq.length); pos++) {
        CalcInside(pos);
        if (pos >= x-1 && pos-1 >= 0) {
            if (ddebug) {
                cout << "in" << pos << " " << Outer(alpha, pos-1) << " " << ReCalcDouterInside(pos) << endl;
                assert(fabs(Outer(alpha, pos-1)-ReCalcDouterInside(pos)) < 0.001);
            }
            Outer(alpha, pos-1) = ReCalcDouterInside(pos);
        }
    }
    for (LEN pos = RightRange(x+1); pos >= max(static_cast<LEN>(0), _start-_constraint*CONST+1); pos--) {
        CalcInsideFromRight(pos);
        if (pos <= x+1 && pos < seq.length) {
            if (ddebug) {
                cout << "out" << pos << " " << Outer(beta, pos) << " " << ReCalcDouterOutside(pos) << endl;
                assert(fabs(Outer(beta, pos)-ReCalcDouterOutside(pos)) < 0.001);
            }
            Outer(beta, pos) = ReCalcDouterOutside(pos);
        }
    }
}

DOUBLE ParasoR::MaxDiff(const Vec& ori, const Vec& mut)
{
    Vec diff;
    transform(ori.begin(), ori.end(), mut.begin(), back_inserter(diff), [&](DOUBLE l, DOUBLE r) {
        return fabs(l-r);
    });
    return *max_element(diff.begin(), diff.end());
}

string ParasoR::GetDiffFile(int mtype, bool acc)
{
    ostringstream os;
    os << STMP << name;
    if (mtype == Arg::Mut::Sub)
        os << "_mut";
    else if (mtype == Arg::Mut::Del)
        os << "_del";
    else
        os << "_ins";
    if (acc) os << "_acc_";
    else os << "_stem_";
    os << id << "_" << chunk << "_" << _constraint;
    return os.str();
}

void ParasoR::CutOffVector(int mtype, Vec& sbppv, LEN pos) {
    sbppv = bppv;
}


void ParasoR::EditVector(int mtype, Vec& ori, Vec& mut, LEN pos)
{
    if (mtype == Arg::Mut::Ins) {
        mut.erase(mut.begin()+pos);
    } else if (mtype == Arg::Mut::Del) {
        ori.erase(ori.begin()+pos-1);
    }
}

void ParasoR::WriteDiffFile(LEN pos, int base, DOUBLE diff, DOUBLE corr, DOUBLE wcorr, string& file, bool app, bool bin)
{
    ofstream ofs;
    if (bin) {
        if (app) ofs.open(file.c_str(), ios::binary|ios::app);
        else ofs.open(file.c_str(), ios::binary|ios::trunc);
        int h = 1;
        if (!app) ofs.write((char*)&h, sizeof(int));
        ofs.write((char*)&pos, sizeof(LEN));
        ofs.write((char*)&base, sizeof(int));
        ofs.write((char*)&diff, sizeof(DOUBLE));
        ofs.write((char*)&corr, sizeof(DOUBLE));
        ofs.write((char*)&wcorr, sizeof(DOUBLE));
    } else {
        if (app) ofs.open(file.c_str(), ios::app);
        else {
            ofs.open(file.c_str(), ios::trunc);
            ofs << "pos\ttype\tmax_diff\tcorrelation\tcorrelation(win)\n";
        }
        ofs << pos << "\t" << base << "\t" << diff << "\t" << corr << "\t" << wcorr << endl;
    }
}

Vec ParasoR::GetSubstrVec(const Vec& vec, LEN pos, int window)
{
    LEN s = max(static_cast<LEN>(0), pos-_start-static_cast<LEN>(window/2));
    LEN e = min(static_cast<LEN>(vec.size()), pos-_start+static_cast<LEN>(window-window/2));
    return Vec(vec.begin()+s, vec.begin()+e);
}

void ParasoR::WriteDiff(int mtype, LEN pos, int base, Vec ori, Vec& mut, bool acc, bool app, bool mout_flag, string& mout, int window)
{
    EditVector(mtype, ori, mut, pos);
    assert(ori.size() == mut.size());
    DOUBLE diff = MaxDiff(ori, mut);
    DOUBLE corr = Correlation(ori, mut);
    DOUBLE wcorr = Correlation(GetSubstrVec(ori, pos, window), GetSubstrVec(mut, pos, window));
    if ((mout_flag && mout == "") || !binary) {
        if (!app)   cout << "pos\ttype\tmax_diff\tcorrelation\tcorrelation(win)\n";
        cout << pos << "\t" << base << "\t" << diff << "\t" << corr << "\t" << wcorr << endl;
    } else if (mout_flag) {
        WriteDiffFile(pos, base, diff, corr, wcorr, mout, app, false);
    } else {
        string file = GetDiffFile(mtype, acc);
        WriteDiffFile(pos, base, diff, corr, wcorr, file, app, true);
    }
}

void ParasoR::ChangeBase(LEN pos, Arg& arg, bool& app)
{
    for (int j = 1; j <= 4; j++) {
        Vec mbppv;
        char temp = (arg.mtype != Arg::Mut::Ins) ? seq.seqget(pos) : 0;
        bool flag = SlidingSequence(pos, j, arg.mtype);
        Recalculation(pos);
        (arg.acc_flag) ? CalcRangeAcc(mbppv, max(0, arg.window), _start-1, _end, arg.mtype)
                       : CalcRangeStem(mbppv, _start-1, _end, arg.mtype);
        if (!flag) continue;
        if (arg.eSDC_flag) {
            ;// WriteeSDC();
        } else {
            WriteDiff(arg.mtype, pos, j, bppv, mbppv, arg.acc_flag, app, arg.mout_flag, arg.mout, arg.window);
            app = true;
        }
        ComeBackChange(pos, arg.mtype, temp);
        if (arg.mtype == Arg::Mut::Del) break;
    }
}
/**
 * Mutates a base from arg.start to arg.end and computes stem probability from _start to _end.
 */
void ParasoR::MutatedStem(Arg& arg)
{
    bool app = false;
    LEN left = max(static_cast<LEN>(0), arg.start-_constraint), right = min(arg.end+_constraint, seq.length);
    _start = left; _end = right; // The range of stem;
    ReadStemVec(bppv, arg.acc_flag);
    if (arg.mtype == Arg::Mut::Ins) _end++;
    if (arg.mtype == Arg::Mut::Del) _end--;
    Vec adouter, bdouter;
    for (LEN i = max(static_cast<LEN>(1), arg.start); i <= arg.end+1; i++) {
        SetOriginalDouter(adouter, bdouter);
        CopyOriginalOuter(adouter, bdouter);
        SlidingMatrixForMut(i, arg.mtype);
        ChangeBase(i, arg, app);
        if (i == arg.end && arg.mtype != Arg::Mut::Ins) break; // Append a character;
    }
    if (!noout && !arg.mout_flag) cout << "#-Wriiten " << GetDiffFile(arg.mtype, arg.acc_flag) << endl;
}

}