#include "part_func.hh"

namespace Rfold {

void ParasoR::SetOriginalDouter(Vec& adouter, Vec& bdouter)
{
    SetRawRangedMatrix(_start);
    adouter = alpha.outer;
    bdouter = beta.outer;
}

void ParasoR::CopyOriginalOuter(Vec& adouter, Vec& bdouter)
{
    alpha.outer = adouter;
    alpha.iend = alpha.istart+alpha.outer.size()-1;
    beta.outer = bdouter;
    beta.iend = beta.istart+beta.outer.size()-1;
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
        alpha.Insert(x);
        beta.Insert(x);
    } else if (mtype == Arg::Mut::Ins) {
        seq.Delete(x);
        alpha.Delete(x);
        beta.Delete(x);
    }
}

bool ParasoR::SlidingSequence(LEN x, int c, int mtype)
{
    if (mtype < 0) return false;
    else if (mtype == Arg::Mut::Sub) {
        return seq.Sub(x, base(c), c);
    } else if (mtype == Arg::Mut::Del) {
        alpha.Delete(x);
        beta.Delete(x);
        return seq.Delete(x);
    } else if (mtype == Arg::Mut::Ins) {
        alpha.Insert(x);
        beta.Insert(x);
        return seq.Insert(x, base(c), c);
    }
    return true;
}


void ParasoR::Recalculation(LEN x)
{
    for (LEN pos = LeftRange(x); pos <= RightRange(_end); pos++) {
        CalcInside(pos);
        if (pos >= x-1 && pos > 0)
            Outer(alpha, pos) = ReCalcDouterInside(pos);
    }
    for (LEN pos = RightRange(x+1); pos >= LeftRange(x); pos--) {
        CalcInsideFromRight(pos);
        if (pos <= x+1 && pos < seq.length)
            Outer(beta, pos) = ReCalcDouterOutside(pos);
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
    os << TMP << name;
    if (mtype == Arg::Mut::Sub)
        os << "_mut_";
    else if (mtype == Arg::Mut::Del)
        os << "_del_";
    else
        os << "_ins_";
    if (acc) os << "_acc_";
    else os << "_stem_";
    os << id << "_" << chunk << "_" << _constraint;
    return os.str();
}

void ParasoR::CutOffVector(int mtype, Vec& sbppv, LEN pos) {}


void ParasoR::EditVector(int mtype, Vec& ori, Vec& mut, LEN pos)
{
    if (mtype == Arg::Mut::Ins) {
        mut.erase(mut.begin()+pos);
    } else if (mtype == Arg::Mut::Del) {
        ori.erase(ori.begin()+pos);
    }
}

void ParasoR::WriteDiffBin(LEN pos, int base, DOUBLE diff, DOUBLE corr, string& file, bool app)
{
    ofstream ofs;
    if (app) ofs.open(file.c_str(), ios::binary|ios::app);
    else ofs.open(file.c_str(), ios::binary|ios::trunc);
    int h = 1;
    if (!app) ofs.write((char*)&h, sizeof(int));
    ofs.write((char*)&pos, sizeof(LEN));
    ofs.write((char*)&base, sizeof(int));
    ofs.write((char*)&diff, sizeof(DOUBLE));
    ofs.write((char*)&corr, sizeof(DOUBLE));
}

void ParasoR::WriteDiff(int mtype, LEN pos, int base, Vec ori, Vec& mut, bool acc, bool app)
{
    assert(ori.size() == mut.size());
    EditVector(mtype, ori, mut, pos);
    DOUBLE diff = MaxDiff(ori, mut);
    DOUBLE corr = Correlation(ori, mut);
    if (binary) {
        string file = GetDiffFile(mtype, acc);
        WriteDiffBin(pos, base, diff, corr, file, app);
    } else {
        cout << pos << "\t" << base << "\t" << diff << "\t" << corr << endl;
    }
}

void ParasoR::ChangeBase(LEN pos, Arg& arg, bool& app)
{
    Vec mbppv, sbppv;
    CutOffVector(arg.mtype, sbppv, pos);
    for (int j = 0; j < 4; j++) {
        char temp = (arg.mtype != Arg::Mut::Ins) ? seq.substr(pos, 1)[0] : '0';
        if (!SlidingSequence(pos, j, arg.mtype)) continue;
        Recalculation(pos);
        (arg.acc_flag) ? CalcRangeAcc(mbppv, max(0, arg.window), pos, arg.mtype) : CalcRangeStem(mbppv, pos, arg.mtype);
        if (arg.eSDC_flag) {
            ;// WriteeSDC();
        } else {
            WriteDiff(arg.mtype, pos, j, bppv, mbppv, arg.acc_flag, app);
            app = true;
        }
        ComeBackChange(pos, arg.mtype, temp);
    }
}

void ParasoR::MutatedStem(Arg& arg)
{
    bool app = false;
    LEN right = arg.end;
    if (arg.mtype == Arg::Mut::Ins && arg.end == seq.length) right++;
    SetRange(LeftRange(arg.start), RightRange(arg.end));
    ReadStemVec(bppv, arg.acc_flag);
    Vec adouter, bdouter;
    SetOriginalDouter(adouter, bdouter);
    for (LEN i = arg.start; i <= right; i++) {
        CopyOriginalOuter(adouter, bdouter);
        SlidingMatrixForMut(i, arg.mtype);
        ChangeBase(i, arg, app);
    }
}

}