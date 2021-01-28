#include "part_func.hh"

namespace Rfold {

DOUBLE ParasoR::GetPairingStruct(LEN i, LEN j, DOUBLE talpha, DOUBLE tbeta)
{
    DOUBLE stem = (j-i == 1) ? 0 : Logsum(Stem(alpha, i, j), SumExtML(seq.slidebp(i+1, j), i, j+1, true, seq));
    DOUBLE value = Logsum(-talpha, stem, -tbeta);
    return value;
}

DOUBLE ParasoR::ReCalcDouterInside(LEN j)
{
    if (j <= 1 || j > seq.length) return 0.0;
    DOUBLE dalpha = 0;
    DOUBLE newDouter = -INF;
    for (LEN k = 1; k <= _constraint+1 && k <= j; k++) {
        DOUBLE tdalpha = Logsum(GetStemDelta(j, j-k, true), -dalpha);
        newDouter = Logsumexp(newDouter, tdalpha);
        if (j-k-1 >= 0) dalpha = Logsum(dalpha, Outer(alpha, j-k-1));
    }
    return newDouter;
}

DOUBLE ParasoR::ReCalcDouterOutside(LEN j)
{
    if (j < 0 || j > seq.length-1) return 0.0;
    DOUBLE dbeta = 0;
    DOUBLE newDouter = -INF;
    for (LEN k = 1; k <= _constraint+1 && j+k <= seq.length; k++) {
        DOUBLE tdbeta = Logsum(GetStemDelta(j, j+k, false), -dbeta);
        newDouter = Logsumexp(newDouter, tdbeta);
        if (j+k <= seq.length) dbeta = Logsum(dbeta, Outer(beta, j+k));
    }
    return newDouter;
}

void ParasoR::PreCalcInside(LEN x)
{
    for (LEN j = max(static_cast<LEN>(0), x-2*_constraint-10); j <= RightRange(x); j++) {
        CalcInside(j);
    }
}

void ParasoR::PreCalcOutside(LEN x)
{
    for (LEN i = LeftRange(x); i < x; i++) {
        CalcOutside(i);
    }
}

void ParasoR::PreCalcOutside(LEN x, DOUBLE uxx)
{
    LEN left = LeftRange(x);
    if (left == 0)
        uxx = 0.0;
    else {
        for (LEN i = x; i > left; i--)
            uxx = RewindLocalOuter(i, uxx);
    }
    // if (ddebug) cout << "#--rewind: " << left << " " << uxx << endl;
    for (LEN i = left; i < x; i++) {
        CalcOutside(i, uxx);
        uxx = SlideLocalOuter(i, uxx);
    }
}

DOUBLE ParasoR::ReproduceLocalOuter(LEN x)
{
    DOUBLE uxx = -INF;
    if (x == 0) return 0.0;
    DOUBLE talpha = 0.0;
    for (LEN i = x; i >= LeftRange(x); i--) {
        DOUBLE tbeta = Outer(beta, x);
        for (LEN j = x+1; j <= RightRange(i); j++) {
            DOUBLE value = GetPairingStruct(i, j, talpha, tbeta);
            uxx = Logsumexp(uxx, value);
            tbeta = tbeta+Outer(beta, j);
        }
        if (i > 0) talpha += Outer(alpha, i-1);
    }
    return -uxx;
}

DOUBLE ParasoR::ExpandLocalOuter(DOUBLE uxx, LEN x, LEN i, LEN j)
{
    for (LEN pos = x-1; pos >= i; pos--) {
        uxx -= Outer(alpha, pos);
    }
    for (LEN pos = x; pos < j; pos++) {
        uxx -= Outer(beta, pos);
    }
    return uxx;
}

DOUBLE ParasoR::ShrinkLocalOuterInside(DOUBLE uij, LEN i) {
    return uij + Outer(alpha, i);
}

DOUBLE ParasoR::ShrinkLocalOuter(DOUBLE uij, LEN j) {
    return uij + Outer(beta, j-1);
}

DOUBLE ParasoR::RewindLocalOuter(LEN pre, DOUBLE old) {
    return old-Outer(alpha, pre-1)+Outer(beta, pre-1);
}

DOUBLE ParasoR::SlideLocalOuter(LEN x, DOUBLE old) {
    // if (ddebug)
    //     cout << "#--slide " << x << " " << Outer(alpha, x) << " " << Outer(beta, x) << endl;
    return old+Outer(alpha, x)-Outer(beta, x);
}

/* ///////////////////////////////////////////// */

void ParasoR::SetOutsideMat(LEN i, LEN j, DOUBLE value)
{
    if (i > 0 && j < seq.length) {
        if (hard || mfe) {
            Stemend(beta, i, j) = GetOutStemendHC(i, j);
            Multi(beta, i, j) = GetOutMultiHC(i, j);
            Multi1(beta, i, j) = GetOutMulti1HC(i, j);
            Multi2(beta, i, j) = GetOutMulti2HC(i, j);
            Multibif(beta, i, j) = GetOutMultiBifHC(i, j);
        } else {
            Stemend(beta, i, j) = GetOutStemend(i, j);
            Multi(beta, i, j) = GetOutMulti(i, j);
            Multi1(beta, i, j) = GetOutMulti1(i, j);
            Multi2(beta, i, j) = GetOutMulti2(i, j);
            Multibif(beta, i, j) = GetOutMultiBif(i, j);
        }
    }
    if (delta) {
        Stem(beta, i, j) = GetOutStem(i, j, value);
    } else if (hard || mfe) {
        Stem(beta, i, j) = GetOutStemHC(i, j);
    } else {
        Stem(beta, i, j) = GetOutStem(i, j);
    }
}

DOUBLE ParasoR::GetOutStem(LEN i, LEN j, DOUBLE value)
{
    LEN p = i+1, q = j, type = seq.slidebp(p, q);
    DOUBLE temp = -INF;
    if (IsPair(type)) {
        temp = Logsum(value, SumExtML(type, i, j+1, true, seq));
        if (q-p-1 >= TURN) {
            for (LEN ip = i; ip >= LeftRange(j); ip--) {
                for (LEN jp = min(ip+_constraint-1, min(j+MAXLOOP-(i-ip), seq.length-1)); jp >= j && (i-ip)+(jp-j) > 0; jp--) {
                    if (IsPair(seq.slidebp(ip, jp+1)))
                        temp = Logsumexp(temp, Logsum(Stemend(beta, ip, jp), LogLoopEnergy(ip, jp+1, p, q, seq)));
                }
            }
        }
        temp = Logsumexp(temp, Logsum(Multi2(beta, i, j), SumExtML(type, i, j+1, false, seq), logMLintern));
        if (IsOnlyRange(i-1, j+1))
            temp = Logsumexp(temp, Logsum(Stem(beta, i-1, j+1), LogLoopEnergy(i, j+1, p, q, seq)));
    }
    return temp;
}


/* ///////////////////////////////////////////// */

int ParasoR::ReadPartConnectedDouter(bool inside, Vec& vec)
{
    string str;
    int tid = (inside) ? ((id == 0) ? id : id-1) : ((id == chunk-1) ? id : id+1);
    string filename = GetShrunkFileList(File::Part, inside, tid);
    if (!noout) cout << "#" << filename << endl;
    ifstream ifs(filename.c_str());
    if (!ifs) return 0;
    if (!noout) cout << "#-Reading " << filename << endl;
    for (int i = 0; inside && i <= tid; i++) {
        if (!getline(ifs, str)) return vec.size();
        if (i == tid) ReadVec(vec, str);
        else {
            Vec temp;
            ReadVec(temp, str);
        }
    }
    for (int i = chunk-1; !inside && i >= tid; i--) {
        if (!getline(ifs, str)) return vec.size();
        if (i == tid) ReadVec(vec, str);
        else {
            Vec temp;
            ReadVec(temp, str);
        }
    }
    return vec.size();
}

int ParasoR::ReadBinPartConnectedDouter(bool inside, Vec& vec)
{
    int h = 0;
    int tid = (inside) ? ((id == 0) ? id : id-1) : ((id == chunk-1) ? id : id+1);
    string filename = GetShrunkFileList(File::Part, inside, tid);
    ifstream ifs(filename.c_str(), ios::binary);
    if (!ifs || (h = GetColumn(ifs)) != 2*(int)_constraint) {
        return 0;
    }
    if (!noout) cout << "#-Reading " << filename << endl;
    if (!noout) cout << "column size: " << h << endl;
    for (LEN i = 0; inside && i <= tid; i++) {
        if (i == tid) ReadBinVec(h, vec, ifs);
        else {
            Vec temp;
            if (!ReadBinVec(h, temp, ifs)) break;
        }
        GetColumn(ifs);
    }
    for (LEN i = chunk-1; !inside && i >= tid; i--) {
        if (i == tid) ReadBinVec(h, vec, ifs);
        else {
            Vec temp;
            if (!ReadBinVec(h, temp, ifs)) break;
        }
        GetColumn(ifs);
    }
    if (!noout) cout << "#--size: " << vec.size() << endl;
    return vec.size();
}

bool ParasoR::ReadConnectedDouter(bool inside)
{
    string str;
    ifstream ifs(GetDoFile(inside).c_str());
    if (!noout) cout << "#-Reading " << GetDoFile(inside) << " " << alpha.istart << "~" << alpha.iend << endl;
    if (!ifs) return false;
    for (LEN i = 0; i <= alpha.iend; i++) {
        if (!getline(ifs, str)) return i == alpha.iend;
        if (str == "") {
            i--; continue;
        }
        if (i >= alpha.istart && str.length() > 0) {
            DOUBLE value = atof(str.c_str());
            if (inside)
                Outer(alpha, i) = value;
            else
                Outer(beta, i) = value;
        }
    }
    return true;
}

bool ParasoR::ReadBinConnectedDouter(bool inside)
{
    ifstream ifs(GetDoFile(inside).c_str(), ios::binary);
    if (!noout) cout << "#-Reading " << GetDoFile(inside) << " " << alpha.istart << "~" << alpha.iend << endl;
    int h = GetColumn(ifs);
    if (!ifs || h != 1) return false;
    DOUBLE value;
    if (!noout) cout << "#--column size: " << h << endl;
    ifs.seekg(sizeof(int)+sizeof(DOUBLE)*(alpha.istart), std::ios::beg);
    for (LEN i = alpha.istart; i <= alpha.iend; i++) {
        ifs.read((char*)&(value), sizeof(DOUBLE));
        if (ifs.eof()) return i == alpha.iend;
        if (inside) Outer(alpha, i) = value;
        else Outer(beta, i) = value;
    }
    return true;
}

void ParasoR::ReadStem(Vec& vec, string file, LEN s, LEN e)
{
    string str;
    ifstream ifs(file.c_str());
    if (!noout) cout << "#-Reading " << file << endl;
    for (LEN i = 0; i < e; i++) {
        if (!getline(ifs, str, '\t')) {
            cerr << file << " does not include enough data\nPlease retry to make " << file << endl;
            exit(-1);
        }
        if (i >= s && str.length() > 0) {
            DOUBLE value = atof(str.c_str());
            vec.push_back(value);
        }
    }
}

void ParasoR::ReadBinStem(Vec& vec, string file, LEN s, LEN e)
{
    ifstream ifs(file.c_str(), ios::binary);
    if (!noout) cout << "#-Reading " << file << endl;
    int h = GetColumn(ifs);
    DOUBLE value;
    if (!noout) cout << "#--column size: " << h << endl;
    if (h != 1) {
        cerr << file << "'s size != 1\nPlease retry to make " << file << endl;
        exit(-1);
    }
    ifs.seekg(sizeof(int)+sizeof(DOUBLE)*(s), std::ios::beg);
    for (LEN i = s; i < e; i++) {
        ifs.read((char*)&(value), sizeof(DOUBLE));
        if (ifs.eof()) break;
        vec.push_back(value);
    }
}

void ParasoR::ReadStemVec(Vec& vec, bool acc)
{
    LEN s = _start;
    LEN e = _end;
    (binary) ? ReadBinStem(vec, GetStemFile(acc), s, e)
             : ReadStem(vec, GetStemFile(acc), s, e);
}

void ParasoR::InitRowMat(Matrix& mat, LEN pos)
{
    InitVec(mat.stem, Ind(mat, pos));
    InitVec(mat.stemend, Ind(mat, pos));
    InitVec(mat.multi, Ind(mat, pos));
    InitVec(mat.multi1, Ind(mat, pos));
    InitVec(mat.multi2, Ind(mat, pos));
    InitVec(mat.multibif, Ind(mat, pos));
}

void ParasoR::InitColMat(Matrix& mat, LEN pos)
{
    for (LEN j = pos; j <= RightRange(pos); j++)
    {
        Stem(mat, pos, j) = -INF;
        Stemend(mat, pos, j) = -INF;
        Multi(mat, pos, j) = -INF;
        Multi1(mat, pos, j) = -INF;
        Multi2(mat, pos, j) = -INF;
        Multibif(mat, pos, j) = -INF;
    }
}

void ParasoR::CalcInside(LEN tpos)
{
    InitRowMat(alpha, tpos);
    for (LEN i = tpos-TURN; i >= LeftRange(tpos); i--) {
        SetInsideMat(i, tpos);
     }
}

void ParasoR::CalcInsideFromRight(LEN tpos)
{
    InitColMat(alpha, tpos);
    for (LEN j = tpos+TURN; j <= RightRange(tpos); j++) {
        SetInsideMat(tpos, j);
    }
}

void ParasoR::CalcOutside(LEN pos)
{
    InitColMat(beta, pos);
    for (LEN j = RightRange(pos); j >= pos+TURN; j--)
        SetOutsideMat(pos, j, 0);
}

void ParasoR::CalcOutside(LEN pos, DOUBLE uxx)
{
    InitColMat(beta, pos);
    LEN right = RightRange(pos);
    DOUBLE uij = ExpandLocalOuter(uxx, pos, pos, right);
    // if (ddebug) cout << "#----expand : " << uxx << endl;
    for (LEN j = right; j >= pos+TURN; j--) {
        // if (ddebug) cout << "#---u" << pos << "_" << j << ": " << uij << endl;
        SetOutsideMat(pos, j, uij);
        uij = ShrinkLocalOuter(uij, j);
    }
}

void ParasoR::CalcForward(LEN pos)
{
    if (seq.length >= pos+_constraint+1)
        CalcInside(pos+_constraint+1);
    CalcOutside(pos);
}

void ParasoR::CalcForward(LEN pos, DOUBLE uxx)
{
    if (seq.length >= pos+_constraint+1)
        CalcInside(pos+_constraint+1);
    CalcOutside(pos, uxx);
}

void ParasoR::SetRawRangedMatrix(bool set)
{
    InitBpp(false, set);
    if (!set) return;
    if (binary) {
        if (!ReadBinConnectedDouter(true) || !ReadBinConnectedDouter(false))
            throw "File read error. Please check connected douter file.";
        if (!noout)
            cout << "#-Finish reading connected douters binary" << endl;
    } else {
        if (!ReadConnectedDouter(true) || !ReadConnectedDouter(false))
            throw "File read error. Please check connected douter file.";
        if (!noout)
            cout << "#-Finish reading connected douters" << endl;
    }
    if (ddebug) {
        cout << "#--inside\n";
        PrintVec(alpha.outer, true);
        cout << "#--outside\n";
        PrintVec(beta.outer, true);
    }
}

DOUBLE ParasoR::SetRangedMatrix(LEN start, bool set)
{
    if (!noout && set) cout << "#-SetMatrix..." << endl;
    SetRawRangedMatrix(set);
    PreCalcInside(min(seq.length-1, start+1));
    DOUBLE uxx = ReproduceLocalOuter(min(seq.length-1, start+1)); // CalcInside(~pos+_constraint);
    PreCalcOutside(min(seq.length-1,start+1), uxx);
    return uxx;
}


void ParasoR::StoreBppSlide(LEN i, LEN right, Mat& bppm)
{
    for (LEN j = max(_start+1, i+TURN+1); j <= i+_constraint && j <= right; j++) {
        // DOUBLE tbpp = (delta) ? bppDelta(i, j, true) : bpp(i, j);
        DOUBLE tbpp = bpp_func(i, j, delta);
        // if (ddebug) cout << "bpp: " << i << " " << j << " " << tbpp << endl;
        BPPM(j, j-i, _start+1) = tbpp;
    }
}

void ParasoR::StoreBppSlide(LEN i, LEN right, Vec& prebpp)
{
    for (LEN j = max(_start+1, i+TURN); j <= i+_constraint && j <= right; j++) {
        // DOUBLE tbpp = (delta) ? bppDelta(i, j, true) : bpp(i, j);
        DOUBLE tbpp = bpp_func(i, j, delta);
        prebpp[i%(_constraint+1)] += tbpp;
        prebpp[j%(_constraint+1)] += tbpp;
    }
}

void ParasoR::StoreAccProfSlide(LEN i, LEN right, int out, bool store)
{
    if (out == Out::ACC) {
        if (i+static_cast<LEN>(_window) <= seq.length) {
            cout << "*\t" << i << "-" << i+static_cast<LEN>(_window) << "\t:\t" << acc(i, i+static_cast<LEN>(_window)) << endl;
        }
    } else if (out == Out::PROF || out == Out::PROFIM) {
        GetProfs(i, prom, 0, 0.0, store);
        if (out == Out::PROF) {
            cout << "*\t" << i << "\t:\t";
            for (LEN j = static_cast<LEN>(prom.size()-TYPE); j < static_cast<LEN>(prom.size()); j++)
                cout << prom[j] << ",";
            cout << endl;
        }
    } else if (out == Out::MOTIF) {
        GetProfs(i, prom);
        auto it = prom.begin();
        cout << substructure(std::distance(it, max_element(it, it+TYPE)));
        if (i == seq.length) cout << endl;
    }
}

void ParasoR::StoreAreaBppSlide(LEN i, LEN right, Mat& bppm)
{
    LEN shift = _start+1;
    for (LEN j = max(_start+1, i+TURN+1); j <= i+_constraint && j <= right; j++) {
        // DOUBLE tbpp = (delta) ? bppDelta(i, j, true) : bpp(i, j);
        DOUBLE tbpp = bpp_func(i, j, delta);
        BPPM(j, j-i, shift) = tbpp;
    }
}

void ParasoR::StoreAreaBppSlide(LEN i, LEN right, Vec& prebpp)
{
    LEN shift = _start+1;
    for (LEN j = max(_start+1, i+TURN+1); j <= i+_constraint && j <= right; j++) {
        // DOUBLE tbpp = (delta) ? bppDelta(i, j, true) : bpp(i, j);
        DOUBLE tbpp = bpp_func(i, j, delta);
        if (i >= shift) prebpp[i-shift] += tbpp;
        if (j <= _end) prebpp[j-shift] += tbpp;
    }
}

void ParasoR::CheckDouter(LEN start, LEN end, DOUBLE uxx)
{
    LEN bstart = start+1;
    LEN right = RightBpRange(end);
    Vec P = Vec(end-start, 0.0);
    for (LEN i = max(static_cast<LEN>(1), bstart); i <= min(seq.length, end+_constraint); ++i)
    {
        DOUBLE diff1 = Outer(alpha, i-1)-ReCalcDouterInside(i);
        DOUBLE diff2 = (i-_constraint > 0) ? Outer(beta, i-_constraint)-ReCalcDouterOutside(i-_constraint) : 0.0;
        if (abs(diff1) > 1e-6 || abs(diff2) > 1e-6) {
            cout << "Diff " << i << " " << diff1 << " " << Outer(alpha, i-1) << endl;
            cout << "Out " << i << " " << diff2 << " " << Outer(beta, i-_constraint) << endl;
        }
        CalcForward(i);
    }
}

void ParasoR::SetProbs(Mat& P) {
    P = Mat(min(seq.length, _end+_constraint+1)-_start, Vec(_constraint+1, -INF));
}
void ParasoR::SetProbs(Vec& P) {
    P = Vec(_end-_start, 0.0);
}

void ParasoR::CalcSlidingWindowBoundary(Vec& P, LEN start, LEN end, bool set)
{
    LEN right = RightBpRange(end);
    LEN bstart = start+1;
    DOUBLE uxx = SetRangedMatrix(start, set);
    // if (ddebug)
    //     CheckDouter(start, end, uxx);
    SetProbs(P);
    if (!noout && set) cout << "#--calculating size:" << right-bstart+1 << endl;
    for (LEN pos = bstart; pos <= end; pos++) {
        P[pos-bstart] = exp(uxx);
        if (!noout && (pos%10000 == 9999)) {
            cout << "#--calculating-- " << pos << endl;
        }
        CalcForward(pos, uxx);
        uxx = SlideLocalOuter(pos, uxx);
    }
}

template <class Probs>
void ParasoR::CalcSlidingWindowStem(Probs& P, LEN start, LEN end, bool set)
{
    LEN right = RightBpRange(end);
    LEN bstart = start+1;
    DOUBLE uxx = SetRangedMatrix(start, set);
    // if (ddebug)
    //     CheckDouter(start, end, uxx);
    SetProbs(P);
    if (!noout && set) cout << "#--calculating size:" << right-bstart+1 << endl;
    for (LEN i = bstart-1; bstart-i <= _constraint && i >= 1; i--) {
        StoreAreaBppSlide(i, right, P);
    }
    for (LEN pos = bstart; pos <= end; pos++) {
        if (!noout && (pos%10000 == 9999)) {
            cout << "#--calculating-- " << pos << endl;
        }
        CalcForward(pos, uxx);
        StoreAreaBppSlide(pos, right, P);
        uxx = SlideLocalOuter(pos, uxx);
    }
}

void ParasoR::CalcSlidingWindowAcc(Vec& P, int region, LEN start, LEN end, bool set)
{
    DOUBLE uxx = SetRangedMatrix(start, set);
    LEN bstart = start+1;
    P = Vec();
    for (LEN pos = bstart; pos <= end; pos++) {
        CalcForward(pos, uxx);
        DOUBLE uij = ExpandLocalOuter(uxx, pos, pos-1, pos+static_cast<LEN>(region));
        // if (ddebug)
        //     cout << "#---au " << pos << " " << pos+static_cast<LEN>(region) << ": " << uij << endl;
        if ((pos-1)%(static_cast<LEN>(region)+1) == 0) {
            P.push_back( (delta) ? accDelta(pos, min(seq.length, pos+static_cast<LEN>(region)), uij)
                                 : acc(pos, min(seq.length, pos+static_cast<LEN>(region))));
        }
        uxx = SlideLocalOuter(pos, uxx);
    }
    if (!noout)
        cout << "#--calculated size " << P.size() << " ("  << start
            << "~" << end << ", " << region+1 << ")" << endl;
}


template <class Structure>
void ParasoR::CalcSlidingWindowProf(Structure& P, LEN start, LEN end, bool set)
{
    LEN right = RightBpRange(end);
    LEN bstart = start+1;
    P = Structure();
    prom = Vec((end-start)*TYPE, 0.);
    DOUBLE uxx = SetRangedMatrix(start, set);
    if (!noout) cout << "#--size:" << right-bstart+1 << endl;
    for (LEN pos = bstart; pos <= end; pos++) {
        if (!noout && (pos%10000 == 0 || pos == bstart)) {
            cout << "#--calculating-- " << pos << endl;
        }
        CalcForward(pos, uxx);
        if (linear) {
            GetProfsLinear(pos, prom, bstart, end, uxx);
        } else {
            GetProfs(pos, prom, bstart, ExpandLocalOuter(uxx, pos, pos-1, pos));
        }
        uxx = SlideLocalOuter(pos, uxx);
    }
    for (LEN pos = end+1; pos <= RightBpRange(end); pos++) {
        if (linear) GetProfsLinearRight(pos, prom, bstart, end);
    }
    AppendProf(P, prom);
}

void ParasoR::CalcStem(Mat& P) {
    CalcSlidingWindowStem(P, _start, _end, true);
}

void ParasoR::CalcStem(Vec& P, bool boundary) {
    if (boundary) CalcSlidingWindowBoundary(P, _start, _end, true);
    else CalcSlidingWindowStem(P, _start, _end, true);
}

void ParasoR::CalcStem(bool store, bool boundary)
{
    Vec P;
    CalcStem(P, boundary);
    if (store && chunk > 0) {
        if (chunk == id)
            StoreStem(GetStemFile(false, false, boundary), P, false);
        else {
            if (!noout)
                cout << "#-Writing to " << GetDividedStemFile(false, false, boundary) << endl;
            StoreStem(GetDividedStemFile(false, false, boundary), P, false);
        }
    } else {
        //cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
        cout.setf(std::ios_base::scientific);
        if (boundary)
            cout << "#--*\tbase\tpos\tboundary prob\n";
        else
            cout << "#--*\tbase\tpos\tstem prob\n";
        for (LEN i = 0; i <= _end-(_start+1); i++)
            cout << "*\t" << setprecision(20) << seq.strget(i+(_start+1)) << "\t" << i+(_start+1) << "\t" << P[i] << endl;
    }
}

void ParasoR::CalcAcc(Vec& P, int region) {
    CalcSlidingWindowAcc(P, region, _start, _end);
}


void ParasoR::CalcAcc(bool store)
{
    Vec P;
    CalcAcc(P, _window-1);
    if (store && chunk > 0) {
        if (chunk == id) {
            StoreStem(GetStemFile(true), P, true);
        } else {
            if (!noout)
                cout << "Writing to " << GetDividedStemFile(true) << endl;
            StoreStem(GetDividedStemFile(true), P, true);
        }
    } else {
        cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
        LEN acc_start = ((_start)/_window + ((_start)%_window > 0 ? 1 : 0))*_window;
        for (LEN i = 0; i < P.size(); i++)
            cout << "*\t" << acc_start+(i*_window)+1 << "~" << acc_start+((i+1)*_window) << "\t" << P[i] << endl;
    }
}

void ParasoR::CalcProf(Vec& P) {
    CalcSlidingWindowProf(P, _start, _end, true);
}
void ParasoR::CalcProf(vector<char>& P) {
    CalcSlidingWindowProf(P, _start, _end, true);
}

void ParasoR::CalcProf(bool value, bool store)
{
    if (value) {
        Vec P;
        CalcProf(P);
        if (store && chunk > 0) {
            if (!noout)
                cout << "#-Writing to " << GetDividedStemFile(true, true) << endl;
            StoreStem(GetDividedStemFile(true, true), P, true);
        } else if (chunk > 0 && chunk == id && store) {
            StoreStem(GetStemFile(true, true), P, true);
        } else {
            cout << "#--pos : Bulge,Outer,Hairpin,Multi,Stem,Interior," << endl;
            for (LEN i = 0; i < P.size(); i++) {
                if (i%TYPE == 0) cout << "\n" << i/TYPE+_start+1 << " : ";
                cout << P[i] << ",";
            }
            cout << endl;
        }
    } else {
        vector<char> P;
        CalcProf(P);
        if (!noout)
            cout << "#-Writing to " << "stdout" << endl;
        for (LEN i = 0; i <= _end-(_start+1) && i < P.size(); i++)
            cout << P[i];
        cout << endl;
    }
}

void ParasoR::GetImage(string str, LEN start, LEN end, Vec stem, int elem)
{
    Plot plot;
    if (end-start <= TURN)  return;
    ostringstream oss;
    oss << STMP << name << "_" << ((elem == PROF_NUM) ? "prof_" : "stem_") << _start << "_se=" << start << "_" << end << "_gamma=" << gamma << ".ps";
    cout << "#-Drawing to " << oss.str() << endl;
    cout << "# local seq: " << seq.substr(start+_start+1, end-start+1) << endl;
    cout << "# local str: " << str.substr(start, end-start+1) << endl;
    Plot::RNAColorPlot(seq.substr(start+_start+1, end-start+1), str.substr(start, end-start+1), oss.str(), stem);
}

bool ParasoR::GetWholeImage(string str, vector<int>& cbpp, Vec& stem, int elem)
{
    Plot plot;
    LEN MAXLEN = 2000;
    if ((_start == static_cast<LEN>(0) && _end == seq.length)) {
        GetImage(str, _start, _end-1, stem, elem);
        return _end-_start < MAXLEN;
    }
    vector<int>::reverse_iterator st = find_if(cbpp.rbegin(), cbpp.rend(), [&str](int x) -> bool {
        return x < 0;
    });
    vector<int>::iterator et = find_if(cbpp.begin(), cbpp.end(), [&str](int x) -> bool {
        return x >= (int)str.length();
    });
    LEN left = (st != cbpp.rend()) ? str.length()-distance(cbpp.rbegin(), st) : 0;
    LEN right = (et != cbpp.end()) ? distance(cbpp.begin(), et)-1 : str.length()-1;
    if (right-left > MAXLEN || right <= left) return false;
    GetImage(str, left, right, Vec(stem.begin()+left*elem, stem.begin()+(right+1)*elem), elem);
    return true;
}

void ParasoR::DrawImage(vector<int>& cbpp, Vec& stem, string& str, int elem)
{
    int left = cbpp.size();
    if (GetWholeImage(str, cbpp, stem, elem)) return;
    for (int i = 0; i < cbpp.size(); i++) {
        if (i > cbpp[i]) {
            left = min(cbpp[i], left);
        } else if (cbpp[i] > i || i-left > _constraint) {
            if (left >= 0 && left < i) {
                GetImage(str, left, cbpp[left], Vec(stem.begin()+left*elem, stem.begin()+(cbpp[left]+1)*elem), elem);
            }
            left = (cbpp[i] > i) ? i : cbpp.size();
            i = cbpp[i];
        }
    }
}

void ParasoR::CalcBpp(bool output, DOUBLE minp, bool calculated)
{
    if (!calculated)
        CalcStem(bppm);
    LEN bstart = _start+1;
    if (!output) return;
    for (LEN j = 0; j < bppm.size(); j++) {
        for (LEN dist = 0; dist < bppm[j].size(); dist++) {
            if (BPPM(j+bstart, dist, bstart) < minp)  continue;
            if (j-dist >= 0 || j < seq.length) cout << "*\t" << j+bstart-dist << "\t" << j+bstart  << "\t" << BPPM(j+bstart, dist, bstart) << endl;
        }
    }
}

void ParasoR::CalcEntropy(bool calculated)
{
    if (!calculated)
        CalcStem(bppm);
    LEN bstart = _start+1;
    for (LEN j = 0; j < min(static_cast<LEN>(bppm.size()), _end-_start); j++) {
        DOUBLE entropy = 0.0;
        DOUBLE stem = 0.0;
        for (LEN dist = 0; dist < bppm[j].size(); dist++) {
            if (BPPM(j+bstart, dist, bstart) > static_cast<DOUBLE>(0.0)) {
                entropy -= ProbEntropy(BPPM(j+bstart, dist, bstart));
                stem += BPPM(j+bstart, dist, bstart);
            }
        }
        for (LEN jp = j+1; jp < bppm.size() && jp+bstart <= RightBpRange(j+bstart); jp++) {
            LEN dist = jp-j;
            if (dist > bppm[jp].size()) continue;
            if (BPPM(jp+bstart, dist, bstart) > static_cast<DOUBLE>(0.0)) {
                entropy -= ProbEntropy(BPPM(jp+bstart, dist, bstart));
                stem += BPPM(jp+bstart, dist, bstart);
            }
        }
        entropy -= ProbEntropy(static_cast<DOUBLE>(1.0-stem));
        cout << "*\t" << seq.strget(j+bstart) << "\t" << j+bstart << "\t" << entropy << endl;
    }
}

void ParasoR::CalcMEA(bool out, bool image, bool prof, bool calculated)
{
    if (!calculated)
        CalcBpp();
    LEN bstart = _start+1;
    bool centroid = (_start == static_cast<LEN>(0) && _end == seq.length);
    if (!noout) {
        if (centroid) cout << "#-Centroid structure (gamma = " << gamma << ")" << endl;
        else cout << "#-Structures having bpp (1/(1+" << gamma << "))" << endl;
    }
    string str =  (centroid) ? Centroid::GetCentroidStructure(bppm, _start, _end, gamma)
    : Centroid::GetMEAStructure(bppm, _end-_start, gamma);
    if (out) cout << "#structure " << str << endl;
    if (image) {
        vector<int> cbpp;
        (centroid) ? Centroid::GetStringBpp(str, cbpp)
        : Centroid::GetMEABpp(bppm, _end-_start, gamma, cbpp);
        if (!prof) {
            Vec stem;
            Centroid::GetStemProb(bppm, _end-_start, stem);
            DrawImage(cbpp, stem, str);
        } else {
            DrawImage(cbpp, prom, str, 6);
        }
    }
}

void ParasoR::CalcRangeStem(Vec& stem, LEN start, LEN end, int mtype)
{
    if (mtype == Arg::Mut::Del) {
        CalcSlidingWindowStem(stem, start, end, false);
    } else if (mtype == Arg::Mut::Ins) {
        CalcSlidingWindowStem(stem, start, end, false);
    } else {
        CalcSlidingWindowStem(stem, start, end, false);
    }
}

void ParasoR::CalcRangeAcc(Vec& stem, int window, LEN start, LEN end, int mtype)
{
    if (mtype == Arg::Mut::Del) {
        CalcSlidingWindowAcc(stem, window, start, end, false);
    } else if (mtype == Arg::Mut::Ins) {
        CalcSlidingWindowAcc(stem, window, start, end, false);
    } else {
        CalcSlidingWindowAcc(stem, window, start, end, false);
    }
}

void ParasoR::TabToNewLine(string& filename)
{
    string oldfile = filename+"_temp";
    string str;
    ifstream ifs(filename.c_str());
    ofstream ofs(oldfile.c_str());
    if (!ifs) return;
    istreambuf_iterator<char> in(ifs), eos;
    ostream_iterator<char> out(ofs);
    copy(in, eos, out);
    ifs.close();
    ofs.close();
    ifs.open(oldfile.c_str());
    ofs.open(filename.c_str());
    for ( ; getline(ifs, str, '\t'); ) {
        ofs << str;
        ofs << "\n";
    }
    ifs.close();
    ofs.close();
    remove(oldfile.c_str());
}

bool ParasoR::ReadStemToSingleFile(string& ifile, string& ofile, bool app, string end)
{
    string str;
    ifstream ifs(ifile.c_str());
    if (!ifs) {
        if (!memory)
            cerr << "No file :" << ifile << endl;
        return true;
    }
    ofstream ofs(ofile.c_str(), ((app) ? ios::app : ios::trunc));
    istreambuf_iterator<char> in(ifs), eos;
    ostream_iterator<char> out(ofs);
    copy(in, eos, out);
    LEN pos = ofs.tellp();
    // if (end != "") {
        ofs.seekp(pos-1);
        ofs << end;
    // }
    ifs.close();
    ofs.close();
    return false;
}

bool ParasoR::ReadBinStemToSingleFile(string& ifile, string& ofile, bool app)
{
    int h = 0;
    ifstream ifs(ifile.c_str(), ios::binary);
    ifs.read((char*)&h, sizeof(int));
    if (!ifs) {
        if (!memory)
            cerr << "No file :" << ifile << endl;
        return true;
    } else if (h < 0) {
        cerr << "File format broken :" << ifile << endl;
        return true;
    } else {
        ofstream ofs;
        if (app) ofs.open(ofile.c_str(), ios::binary|ios::app);
        else ofs.open(ofile.c_str(), ios::binary|ios::trunc);
        h = 1;
        if (!app) ofs.write((char*)&h, sizeof(int));
        istreambuf_iterator<char> in(ifs), eos;
        ostream_iterator<char> out(ofs);
        copy(in, eos, out);
        ifs.close();
        ofs.close();
        return false;
    }
}

void ParasoR::ConcatStemdb(bool acc, bool prof)
{
    if (prof && !acc) return;
    string outfile = GetStemFile(acc, prof);
    for (int i = 0; i < chunk; i++) {
        id = i;
        string file = GetDividedStemFile(acc, prof);
        cout << "#-Reading " << file << endl;
        bool error = (binary) ? ReadBinStemToSingleFile(file, outfile, i != 0)
                 : ReadStemToSingleFile(file, outfile, i != 0, (i == chunk-1) ? "\n" : "\t");
        if (error) return;
    }
    if (!noout) cout << "#-Concatenation succeeded. (-> " << outfile << ")" << endl;
    RemoveStem(acc, prof);
}

void ParasoR::Connect(ParasoR& rfold, bool shrink, bool keep)
{
    if (rfold.memory) {
        if (!noout) cout << "#--Memory saving mode" << endl;
        if (!rfold.ConnectSavedFiles(keep))
            rfold.ConnectDoSaved(keep);
    } else {
        if (!noout) cout << "#--file: " << rfold.GetDoFile(true) << ", " << rfold.GetDoFile(false) << endl;
        rfold.ConnectDo(keep);
        rfold.ConcatDo();
    }
}

void ParasoR::Connect(Arg& arg, bool shrink)
{
    if (arg.str.length() == 0 && arg.length == 0) return;
    ParasoR rfold;
    rfold.SetBasicParam(arg, shrink);
    if (!rfold.SetChunkId(0, arg.chunk)) return;
    if (arg.init_file) {
        rfold.RemoveTempFilesSaveMem(true);
        rfold.RemoveTempFilesSaveMem(false);
        cout << "Complete file initiation." << endl;
        return;
    }
    Connect(rfold, shrink, arg.keep_flag);
    if (arg.stemdb_flag) {
        rfold.SetRange(1, arg.length);
        (arg.acc_flag) ? rfold.CalcAcc(true) : rfold.CalcStem(true);
    }
}

void ParasoR::Stemdb(Arg& arg, bool shrink)
{
    if (arg.str.length() == 0 && arg.length == 0) return;
    ParasoR rfold;
    rfold.SetBasicParam(arg, shrink);
    rfold.cut = true;
    rfold.SetWindow(arg.window);
    if (arg.id == arg.chunk) {
        rfold.SetChunkId(0, arg.chunk);
        rfold.ConcatStemdb(arg.acc_flag, arg.prof_flag);
    } else if (rfold.SetChunkId(arg.id, arg.chunk)) {
        if (arg.acc_flag && !arg.prof_flag) {
            rfold.CalcAcc(!arg.outtext);
        } else if (arg.mea_flag) {
            rfold.CalcMEA(true, arg.image, arg.prof_flag);
        } else if (arg.prof_flag) {
            rfold.CalcProf(arg.acc_flag, !arg.outtext);
        } else if (arg.entro_flag) {
            rfold.CalcEntropy(!arg.outtext);
        } else {
            rfold.CalcStem(!arg.outtext, arg.boundary);
        }
    } else {
        cerr << "too many chunk for this sequence len: " << rfold.seq.length << endl;
    }
}

/* ///////////////////////////////////////////// */

void ParasoR::CalcInsideOuter(LEN j)
{
    assert(delta == false);
    if (j == 0) Outer(alpha, 0) = 0.0;
    else {
        Outer(alpha, j) = Outer(alpha, j-1);
        DOUBLE temp = -INF;
        for (LEN k = LeftRange(j); k < j; k++) {
            temp = Logsumexp(temp, Logsum(Outer(alpha, k), Stem(alpha, k, j), SumExtML(seq.slidebp(k+1, j), k, j+1, true, seq)));
            // if (debug && !Is_INF(Stem(alpha, k, j))) {
            //     cout << "outer: " << j << endl;
            //     cout << Outer(alpha, k) << " " << Stem(alpha, k, j) << " " << SumExtML(seq.slidebp(k+1, j), k, j+1, true, seq) << endl;
            // }
        }
        Outer(alpha, j) = Logsumexp(temp, Outer(alpha, j));
    }
}

void ParasoR::CalcOutsideOuter(LEN j)
{
    assert(delta == false);
    if (j == seq.length) Outer(beta, seq.length) = 0.0;
    else {
        Outer(beta, j) = Outer(beta, j+1);
        DOUBLE temp = -INF;
        for (LEN k = j+TURN; k <= RightRange(j); k++) {
            temp = Logsumexp(temp, Logsum(Stem(alpha, j, k), Outer(beta, k), SumExtML(seq.slidebp(j+1, k), j, k+1, true, seq)));
        }
        Outer(beta, j) = Logsumexp(Outer(beta, j), temp);
    }
}


void ParasoR::CalcOuter()
{
    for (LEN j = TURN-1; j <= seq.length; j++) {
        CalcInside(j);
        CalcInsideOuter(j);
    }
    for (LEN i = seq.length-TURN; i >= 0; i--) {
        CalcInsideFromRight(i);
        CalcOutsideOuter(i);
    }
}

void ParasoR::CalcBppAtOnce(int out, bool image, DOUBLE thres)
{
    SetProbs(bppm);
    for (LEN pos = 1, right = RightBpRange(seq.length); pos <= seq.length; pos++) {
        CalcForward(pos);
        StoreBppSlide(pos, seq.length, bppm);
    }
    if (out == Out::BPPIM)
        CalcMEA(true, image, false, true);
    else if (out == Out::PROFIM)
        CalcMEA(true, image, true, true);
    else if (out == Out::BPP)
        CalcBpp(true, thres, true);
    else if (out == Out::ENTRO)
        CalcEntropy(true);
}

void ParasoR::PrintRfoldHeader(int out)
{
    cout << std::setprecision(10);
    cout << "#--log(Z): " << alpha.outer[seq.length] << endl;
    if (mfe)
        cout << "#--mfe-enabled" << endl;;
    if (hard)
        cout << "#--hard-const-enabled" << endl;
    else
        cout << "#--hard-const-disabled" << endl;
    if (out == Out::ACC)
        cout << "#--range : accessibility (kcal/mol)" << endl;
    else if (out == Out::PROF)
        cout << "#--pos : Bulge,Outer,Hairpin,Multi,Stem,Interior," << endl;
    else if (out == Out::MOTIF) {
        cout << "#-- B(ulge) or O(uter) or H(airpin) or M(ulti) or S(tem) or I(nterior)" << endl;
        cout << "> " << name << endl;
    } else if (out == Out::ENTRO) {
        cout << "#--entropy sum -(bpp)log2(bpp) " << endl;
    }
}

void ParasoR::CalcAllAtOnce(int out, bool image, DOUBLE thres, bool store)
{
    delta = false;
    InitBpp(false);
    (hard || mfe) ? CalcOuterHC() : CalcOuter();
    Vec prebpp = Vec(_constraint+1, 0);
    _start = 0;
    PrintRfoldHeader(out);
    PreCalcInside(1);
    PreCalcOutside(1);
    if (out == Out::BPP || out == Out::BPPIM || out == Out::ENTRO)
        CalcBppAtOnce(out, image, thres);
    else {
        for (LEN pos = 1, right = RightBpRange(seq.length); pos <= seq.length; pos++) {
            CalcForward(pos);
            if (out == Out::STEM) {
                StoreBppSlide(pos, seq.length, prebpp);
                if (pos%JUMP == 0) {
                    cout << "*\t" << seq.strget(pos) << "\t:\t" << prebpp[pos%(_constraint+1)] << endl;
                }
                prebpp[pos%(_constraint+1)] = 0.0;
            } else {
                StoreAccProfSlide(pos, seq.length, out, out == Out::PROFIM);
            }
        }
    }
    if (out == Out::PROFIM) {
        PreCalcInside(1);
        PreCalcOutside(1);
        CalcBppAtOnce(out, image, thres);
    }
}


}
