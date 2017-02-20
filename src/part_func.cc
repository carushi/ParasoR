#include "part_func.hh"

namespace Rfold {

DOUBLE ParasoR::GetInStem(LEN i, LEN j)
{
    int bp1 = seq.slidebp(i+1, j);
    DOUBLE temp = -INF;
    if (IsPair(bp1)) {
        if (!Is_INF(Stem(alpha, i+1, j-1))) {
            temp = Logsumexp(temp, Logsum(Stem(alpha, i+1, j-1), LogLoopEnergy(i+1, j, i+2, j-1, seq)));
        }
        temp = Logsumexp(temp, Stemend(alpha, i+1, j-1));
    }
    return temp;
}

DOUBLE ParasoR::GetInMultiBif(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    for (LEN k = i+TURN+1; k < j-TURN; k++)
        temp = Logsumexp(temp, Logsum(Multi1(alpha, i, k), Multi2(alpha, k, j)));
    return temp;
}

DOUBLE ParasoR::GetInMulti2(LEN i, LEN j)
{
    DOUBLE temp = Logsum(Multi2(alpha, i, j-1), logML_BASE);
    int type = seq.slidebp(i+1, j);
    if (IsPair(type)) {
        temp = Logsumexp(temp, Logsum(Stem(alpha, i, j), SumExtML(type, i, j+1, false, seq), logMLintern));
    }
    return temp;
}

DOUBLE ParasoR::GetInMulti1(LEN i, LEN j) {
    DOUBLE temp = Logsumexp(Multi2(alpha, i, j), Multibif(alpha, i, j));
    return temp;
}

DOUBLE ParasoR::GetInMulti(LEN i, LEN j) {
    DOUBLE temp = Logsumexp(Multibif(alpha, i, j),
                            Logsum(Multi(alpha, i+1, j), logML_BASE));
    return temp;
}

DOUBLE ParasoR::GetInStemend(LEN i, LEN j)
{
    LEN p = i, q = j+1;
    DOUBLE temp = -INF;
    if (q > p+_constraint) return temp;
    else if (q > seq.length || !IsPair(seq.slidebp(p, q))) return temp;
    else if (p > 0 && q <= seq.length) {
        temp = LogHairpinEnergy(p, q, seq);
        for (LEN ip = i; ip < j-TURN-1; ip++) {
            for (LEN jp = max(j-MAXLOOP+(ip-i), ip+TURN+2); jp <= j && (j-jp)+(ip-i) > 0; jp++) {
                if (IsPair(seq.slidebp(ip+1, jp))) {
                    temp = Logsumexp(temp, Logsum(Stem(alpha, ip, jp), LogLoopEnergy(p, q, ip+1, jp, seq)));
                }
            }
        }
        temp = Logsumexp(temp, Logsum(Multi(alpha, i, j), SumExtML(seq.sliderbp(p, q), j, i+1, false, seq),
                                      logMLclosing+logMLintern));
    }
    return temp;
}

void ParasoR::SetInsideMat(LEN i, LEN j)
{
    Stem(alpha, i, j) = GetInStem(i, j);
    Multibif(alpha, i, j) = GetInMultiBif(i, j);
    Multi2(alpha, i, j) = GetInMulti2(i, j);
    Multi1(alpha, i, j) = GetInMulti1(i, j);
    Multi(alpha, i, j) = GetInMulti(i, j);
    Stemend(alpha, i, j) = GetInStemend(i, j);
}

/* ///////////////////////////////////////////// */

DOUBLE ParasoR::GetOutStemend(LEN i, LEN j) {
    return (IsOnlyRange(i-1, j+1) ? Stem(beta, i-1, j+1) : -INF);
}

DOUBLE ParasoR::GetOutMulti(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    if (IsOnlyRange(i-1, j))
        temp = Logsum(Multi(beta, i-1, j), logML_BASE);
    if (IsOnlyRange(i, j+1))
        temp = Logsumexp(temp, Logsum(Stemend(beta, i, j),
                                      SumExtML(seq.sliderbp(i, j+1), j, i+1, false, seq),
                                      logMLclosing+logMLintern));
    return temp;
}

DOUBLE ParasoR::GetOutMulti1(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    for (LEN k = j+TURN+2; k < min(seq.length, i+_constraint+2); k++)
        temp = Logsumexp(temp, Logsum(Multibif(beta, i, k), Multi2(alpha, j, k)));
    return temp;
}

DOUBLE ParasoR::GetOutMulti2(LEN i, LEN j)
{
    DOUBLE temp = Multi1(beta, i, j);
    if (IsOnlyRange(i, j+1))
        temp = Logsumexp(temp, Logsum(Multi2(beta, i, j+1), logML_BASE));
    for (LEN k = max(static_cast<LEN>(0), j-_constraint-1); k < i-TURN-1; k++)
        temp = Logsumexp(temp, Logsum(Multibif(beta, k, j), Multi1(alpha, k, i)));
    return temp;
}

DOUBLE ParasoR::GetOutMultiBif(LEN i, LEN j) {
    return Logsumexp(Multi1(beta, i, j), Multi(beta, i, j));
}

DOUBLE ParasoR::GetOutStem(LEN i, LEN j)
{
    return GetOutStem(i, j, Logsum(Outer(alpha, i), Outer(beta, j)));
    // LEN p = i+1, q = j, type = seq.slidebp(p, q);
    // DOUBLE temp = -INF;
    // if (IsPair(type)) {
    //     temp = Logsum(Outer(alpha, i), Outer(beta, j), SumExtML(type, i, j+1, true, seq));
    //     if (q-p-1 >= TURN) {
    //         for (LEN ip = i; ip >= LeftRange(j); ip--) {
    //             for (LEN jp = min(ip+_constraint-1, min(j+MAXLOOP-(i-ip), seq.length-1)); jp >= j && (i-ip)+(jp-j) > 0; jp--) {
    //                 if (IsPair(seq.slidebp(ip, jp+1))) {
    //                     temp = Logsumexp(temp, Logsum(Stemend(beta, ip, jp), LogLoopEnergy(ip, jp+1, p, q, seq)));
    //                 }
    //             }
    //         }
    //     }
    //     temp = Logsumexp(temp, Logsum(Multi2(beta, i, j), SumExtML(type, i, j+1, false, seq), logMLintern));
    //     if (IsOnlyRange(i-1, j+1)) {
    //         temp = Logsumexp(temp, Logsum(Stem(beta, i-1, j+1), LogLoopEnergy(i, j+1, p, q, seq)));
    //     }
    // }
    // return temp;
}

void ParasoR::SetOutsideMat(LEN i, LEN j)
{
    if (i > 0 && j < seq.length) {
        Stemend(beta, i, j) = GetOutStemend(i, j);
        Multi(beta, i, j) = GetOutMulti(i, j);
        Multi1(beta, i, j) = GetOutMulti1(i, j);
        Multi2(beta, i, j) = GetOutMulti2(i, j);
        Multibif(beta, i, j) = GetOutMultiBif(i, j);
    }
    Stem(beta, i, j) = GetOutStem(i, j);
}

/* ///////////////////////////////////////////// */

DOUBLE ParasoR::GetStemDelta(LEN j, LEN i, bool inside)
{
    if (inside) {
        if (j-i == 1) return 0;     // extension of outer;
        else return Logsum(Stem(alpha, i, j), SumExtML(seq.slidebp(i+1, j), i, j+1, true, seq));
    } else {
        if (i-j == 1) return 0;     // extension of outer;
        else return Logsum(Stem(alpha, j, i), SumExtML(seq.slidebp(j+1, i), j, i+1, true, seq));
    }
}

DOUBLE ParasoR::GetStemDelta(LEN j, LEN i, int h, bool inside)
{
    if (inside) {
        if (i >= _start+1 || i == _start-static_cast<LEN>(h)) {
            if (j-i == 1) return 0;     // extension of outer;
            else return Logsum(Stem(alpha, i, j), SumExtML(seq.slidebp(i+1, j), i, j+1, true, seq));
        }
    } else {
        if (i <= _end-1 || i == _end+static_cast<LEN>(h)) {
            if (i-j == 1) return 0;     // extension of outer;
            else return Logsum(Stem(alpha, j, i), SumExtML(seq.slidebp(j+1, i), j, i+1, true, seq));
        }
    }
    return -INF;
}


DOUBLE ParasoR::CalcDalpha(LEN j, LEN i, int h, DOUBLE dalpha)
{
    DOUBLE value = GetStemDelta(j, i, h, true);
    if (Is_INF(value)) return -INF;
    // if (ddebug) {
    //     cout << value << " " << Stem(alpha, i, j) << " " << SumExtML(seq.slidebp(i+1, j), i, j+1, true, seq)
    //     << " " << GetDo(alpha, i, h) << " " << dalpha << endl;
    //     if (j-i != 1) cout << seq.strget(i+1) << " " << seq.strget(j) << endl;
    // }
    return Logsum(value, GetDo(alpha, i, h), -dalpha);
}

DOUBLE ParasoR::CalcDbeta(LEN j, LEN i, int h, DOUBLE dbeta)
{
    DOUBLE value = GetStemDelta(j, i, h, false);
    if (Is_INF(value)) return -INF;
    // if (ddebug) {
    //     cout << value << " " << Stem(alpha, j, i) << " " << SumExtML(seq.slidebp(j+1, i), j, i+1, true, seq)
    //     << " " << GetDo(beta, i, h) << " " << dbeta << endl;
    //     if (i-j != 1) cout << seq.strget(j+1) << " " << seq.strget(i) << endl;
    // }
    return Logsum(Logsum(value, GetDo(beta, i, h)), -dbeta);
}


void ParasoR::CalcInDeltaOuter(LEN j)
{
    cout.precision(15);
    // if (ddebug) cout << j << endl;
    for (LEN h = 0; h <= _constraint; h++) {
        // if (ddebug) cout << "#------------" << h << endl;
        DOUBLE dalpha = 0, temp = -INF;
        for (LEN k = 1; k <= _constraint+1 && j-k >= 0; k++) {
            dalpha = Logsum(dalpha, GetDo(alpha, j-k, 0));
            if (j-k == _start-h || j >= _start+1) {
                DOUBLE tdalpha = CalcDalpha(j, j-k, h, dalpha);
                temp = Logsumexp(temp, tdalpha);
                // if (ddebug) {
                //     if (j-k == _start-h)
                //         cout << "start-h " << tdalpha << " " << j-k << "=" << _start << "-" << h << endl;
                //     cout << "value " << tdalpha << " " << GetDo(alpha, j-k, h) << " " << Stem(alpha, j-k, j) << endl;
                // }
            }
            // if (ddebug && !Is_INF(temp)) cout << k << ": " << temp << " " << dalpha << endl;
        }
        DOuter(alpha, j, h) = temp;
    }
    // if (ddebug) cout << endl;
}

void ParasoR::CalcChunkInside(bool outer)
{
    // if (ddebug) cout << "inside" << endl;
    for (LEN j = max(_start-_constraint, static_cast<LEN>(0)); j <= _end; j++) {
        // if (ddebug) cout << j << endl;
        CalcInside(j);
        if (outer && j >= _start+1) CalcInDeltaOuter(j);
    }
    // if (ddebug) PrintMat(true);
}

void ParasoR::CalcOutDeltaOuter(LEN j)
{
     for (LEN h = 0; h < _constraint+1; h++) {
        // if (ddebug) cout << "#------------" << h << endl;
        DOUBLE dbeta = 0;
        DOUBLE temp = -INF;
        for (LEN k = 1; j+k <= RightRange(j); k++) {
            dbeta = Logsum(dbeta, GetDo(beta, j+k, 0));
            if (j+k == _end+h || j+k <= _end-1) {
                DOUBLE tdbeta = CalcDbeta(j, j+k, h, dbeta);
                temp = Logsumexp(temp, tdbeta);
                // if (ddebug) {
                //     if (j+k == _end+h)
                //         cout << "end+h+1 " << tdbeta << " " << j+k << "=" << _end << "+" << h << "+" << endl;
                //     cout << "value " << tdbeta << " " << GetDo(alpha, j+k, h) << " " << Stem(alpha, j, j+k) << endl;
                // }
            }
            // if (ddebug && !Is_INF(temp)) cout << k << ": " << temp << " " << dbeta << endl;
        }
        DOuter(beta, j, h) = temp;
    }
    // if (ddebug) cout << endl;
}

void ParasoR::CalcChunkOutside()
{
    // if (ddebug) cout << "outside" << endl;
    for (LEN i = RightRange(_end); i >= _start; i--) {
        // if (ddebug) cout << "i: " << i << endl;
        CalcInsideFromRight(i);
        if (i <= _end-1) CalcOutDeltaOuter(i);
    }
    // if (ddebug) PrintMat(true);
}

void ParasoR::CalcDeltaInOut(bool inside)
{
    if (inside) CalcChunkInside();
    else CalcChunkOutside();
    if (memory) {
        Vec old_douter = Vec();
        int icount = (binary) ? ReadBinPartConnectedDouter(inside, old_douter) : ReadPartConnectedDouter(inside, old_douter);
        if (icount > 0) {
            ConnectDoSaved(inside, old_douter);
        } else {
            StoreDouterTempShrunk(inside);
        }
    } else {
        StoreDouterTemp(inside);
    }
}

void ParasoR::CalcDeltaInOut()
{
    CalcDeltaInOut(true);
    CalcDeltaInOut(false);
}

/* ///////////////////////////////////////////// */

void ParasoR::ReadVec(Vec& tdouter, string& str)
{
    istringstream iss(str);
    vector<string> words;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(words));
    for (vector<string>::iterator it = words.begin(); it != words.end(); it++) {
        tdouter.push_back(atof(it->c_str()));
    }
}

bool ParasoR::ReadBinVec(int h, Vec& tdouter, ifstream& ifs)
{
    DOUBLE value;
    for (int i = 0; i < h; i++) {
        ifs.read((char*)&value, sizeof(DOUBLE));
        if (ifs.eof()) return false;
        tdouter.push_back(value);
    }
    return true;
}

bool ParasoR::ReadDouterInside(Mat& douter, string filename, Vec& first_douter)
{
    string str;
    ifstream ifs(filename.c_str());
    if (!ifs) return false;
    if (!noout) cout << "#-Reading " << filename << endl;
    for (LEN i = 0; getline(ifs, str); i++) {
        if (first_douter.size() > 0 && i == 0) {
            douter.push_back((first_douter));
        } else {
            douter.push_back(Vec());
            ReadVec(douter[i], str);
        }
    }
    if (!noout) cout << "#--size: " << douter.size() << endl;
    return douter.size() > 0;
}

bool ParasoR::ReadBinDouterInside(Mat& douter, string filename, Vec& first_douter)
{
    int h = 0;
    ifstream ifs(filename.c_str(), ios::binary);
    if (!ifs || (h = GetColumn(ifs)) <= 0) return false;
    if (!noout) cout << "#-Reading " << filename << endl;
    if (!noout) cout << "column size: " << h << endl;
    for (LEN i = 0; ; i++) {
        Vec temp = Vec();
        if (ReadBinVec(h, temp, ifs)) {
            if (first_douter.size() > 0 && i == 0) {
                douter.push_back(first_douter);
            } else douter.push_back(temp);
        } else  break;
    }
    if (!noout) cout << "#--size: " << douter.size() << endl;
    return douter.size() > 0;
}

bool ParasoR::ReadDouterOutside(Mat& douter, string filename, Vec& first_douter)
{
    string str;
    ifstream ifs(filename.c_str());
    if (!noout) cout << "#-Reading " << filename << endl;
    if (!ifs) return false;
    for (LEN i = 0; getline(ifs, str); i++) {
        douter.push_back(Vec());
        ReadVec(douter[i], str);
    }
    if (first_douter.size() > 0 && douter.size() > 0)
        douter[static_cast<LEN>(douter.size()-1)] = first_douter;
    if (!noout) cout << "#--size: " << douter.size() << endl;
    return douter.size() > 0;
}

bool ParasoR::ReadBinDouterOutside(Mat& douter, string filename, Vec& first_douter)
{
    int h = 0;
    ifstream ifs(filename.c_str(), ios::binary);
    if (!ifs || (h = GetColumn(ifs)) <= 0) return false;
    if (!noout) cout << "#-Reading " << filename << endl;
    if (!noout) cout << "column size: " << h << endl;
    for ( ; ; ) {
        Vec temp = Vec();
        if (ReadBinVec(h, temp, ifs)) {
            douter.push_back(temp);
        } else  break;
    }
    if (first_douter.size() > 0 && douter.size() > 0)
        douter[static_cast<LEN>(douter.size()-1)] = first_douter;
    if (!noout) cout << "#--size: " << douter.size() << endl;
    return douter.size() > 0;
}


void ParasoR::GetSumDouterList(const Vec& old_douter, Vec& sum_douter, bool inside)
{
    if (inside) {
        DOUBLE acc = 0.0;
        for (LEN i = 0; i <= _constraint; i++) {
            sum_douter[i] = acc;
            if (i < _constraint && i < static_cast<LEN>(old_douter.size()))
                acc = Logsum(acc, old_douter[static_cast<LEN>(old_douter.size()-1)-i]);
        }
    } else {
        DOUBLE acc = 0.0;
        for (LEN i = 0; i <= _constraint; i++) {
            sum_douter[i] = acc;
            if (i < _constraint && i < (int)old_douter.size())
                acc = Logsum(acc, old_douter[i]);
        }
    }
}

DOUBLE ParasoR::GetDenominator(const Vec& douter, const Vec& sum_douter, bool inside)
{
    if (inside) {
        DOUBLE value = -INF;
        for (LEN h = 0; h <= _constraint && h < static_cast<LEN>(sum_douter.size()); h++) {
            if (!Is_INF(sum_douter[h]))
            value = Logsumexp(value, douter[h]-sum_douter[h]);
            // if (ddebug) cout << "deno " << douter[h] << " " << sum_douter[h] << endl;
        }
        return value;
    } else {
        DOUBLE value = -INF;
        for (LEN h = 0; h <= _constraint && h < static_cast<LEN>(sum_douter.size()); h++) {
            if (!Is_INF(sum_douter[h]))
            value = Logsumexp(value, douter[h]-sum_douter[h]);
            // if (ddebug) cout << "deno " << douter[h] << " " << sum_douter[h] << endl;
        }
        return value;
    }
}

DOUBLE ParasoR::GetNumerator(const Mat& douter, const Vec& sum_douter, LEN k, bool inside)
{
    if (inside) {
        DOUBLE value = -INF;
        for (LEN h = 0; h <= _constraint && h < static_cast<LEN>(sum_douter.size()); h++) {
            // if (ddebug) cout << "*****" << douter[k][h] << " " << douter[k-1][0] << " " << sum_douter[h] << endl;
            value = Logsumexp(value, douter[k][h]+douter[k-1][0]-sum_douter[h]);
        }
        return value;
    } else {
        DOUBLE value = -INF;
        for (LEN h = 0; h <= _constraint && h < static_cast<LEN>(sum_douter.size()); h++) {
            // if (ddebug) cout << "*****" << douter[k][h] << " " << douter[k+1][0] << " " << sum_douter[h] << endl;
            value = Logsumexp(value, douter[k][h]+douter[k+1][0]-sum_douter[h]);
        }
        return value;
    }
}

void ParasoR::ConnectInDo(Vec& old_douter, Mat& douter, int tid, string filename, bool app, bool start)
{
    Vec new_douter = Vec(douter.size(), 0.0);
    Vec sum_douter = Vec(old_douter.size(), 0.0);
    if (tid != 0)
        GetSumDouterList(old_douter, sum_douter, true);
    for (LEN i = 1; i < new_douter.size(); i++) {
        // if (ddebug) cout << "#--------" << i << " " << seq.strget(seq.length/chunk*tid+i)<< endl;
        DOUBLE denominator = (i == 1 && start) ? douter[i-1][0] : //start position;
                                      GetDenominator(douter[i-1], sum_douter, true);
        DOUBLE numerator = GetNumerator(douter, sum_douter, i, true);
        new_douter[i] = numerator-denominator;
        // if (ddebug)
        //     cout << new_douter[i] << " " << numerator << " " << denominator << endl;
    }
    StoreDouter(filename, new_douter, true, app);
    old_douter = new_douter;
}

void ParasoR::ConnectOutDo(Vec& old_douter, Mat& douter, int tid, string filename, bool app, bool start)
{
    Vec new_douter = Vec(douter.size(), 0.0);
    Vec sum_douter = Vec(old_douter.size(), 0.0);
    if (tid != chunk-1)
        GetSumDouterList(old_douter, sum_douter, false);
    for (LEN i = douter.size()-2; i >= 0; i--) {
        // if (ddebug) cout << "#--------" << i << " " << seq.strget(seq.length/chunk*tid+i+1) << endl;
        DOUBLE denominator = (i == douter.size()-2 && start) ? douter[i+1][0] : // start position;
                             GetDenominator(douter[i+1], sum_douter, false);
        DOUBLE numerator = GetNumerator(douter, sum_douter, i, false);
        new_douter[i] = numerator-denominator;
        // if (ddebug)
        //     cout << new_douter[i] << " " << numerator << " " << denominator << endl;
    }
    StoreDouter(filename, new_douter, false, app);
    old_douter = new_douter;
}

bool ParasoR::ConnectDo(bool keep_flag, bool inside)
{
    bool flag = true;
    Init(false, true);
    Vec old_douter = Vec(1, 0);
    Vec first_douter;
    for (int i = 0; i < chunk; i++) {
        int tid = (inside) ? i : chunk-1-i;
        Mat douter;
        string file = GetTempFileList(inside, tid);
        if (inside) {
            flag = (binary) ? ReadBinDouterInside(douter, file, first_douter)
                    : ReadDouterInside(douter, file, first_douter);
        } else {
            flag = (binary) ? ReadBinDouterOutside(douter, file, first_douter)
                    : ReadDouterOutside(douter, file, first_douter);
        }
        if (!flag) {
            if (!noout) cout << "File format error" << endl;
            return false;
        }
        if (inside) {
            ConnectInDo(old_douter, douter, tid, GetDoFile(inside), i != 0);
            first_douter = douter[static_cast<LEN>(douter.size()-1)];
        } else {
            ConnectOutDo(old_douter, douter, tid, GetDoFile(inside), i != 0);
            first_douter = douter[0];
        }
    }
    if (!keep_flag) RemoveTempFiles(inside);
    if (!noout) cout << "#-Complete " << GetDoFile(inside) << endl;
    return true;
}

void ParasoR::ConnectDo(bool keep_flag)
{
    Init(false, true);
    ConnectDo(keep_flag, true);
    ConnectDo(keep_flag, false);
}

void ParasoR::ConnectDoSaved(bool keep_flag)
{
    bool iflag = true, oflag = true;
    Init(false, true);
    Vec old_douter = Vec(1, 0);
    Vec first_douter;
    for (int i = 0; i < chunk; i++) {
        Mat douter;
        string ifile = GetShrunkFileList(File::Shrunk, true, i);
        iflag = (binary) ? ReadBinDouterInside(douter, ifile, first_douter)
                : ReadDouterInside(douter, ifile, first_douter);
        if (iflag) {
            ConnectInDo(old_douter, douter, i, GetShrunkFileList(File::Part, true, i), i != 0, false);
            first_douter = douter[static_cast<LEN>(douter.size())-1];
        } else break;
    }
    old_douter = Vec(1, 0);
    first_douter.clear();
    for (int i = chunk-1; i >= 0; i--) {
        Mat douter;
        string ofile = GetShrunkFileList(File::Shrunk, false, i);
        oflag = (binary) ? ReadBinDouterOutside(douter, ofile, first_douter)
                : ReadDouterOutside(douter, ofile, first_douter);
        if (oflag) {
            ConnectOutDo(old_douter, douter, i, GetShrunkFileList(File::Part, false, i), i != chunk-1, false);
            first_douter = douter[0];
        } else break;
    }
}

void ParasoR::ReadBinFirstDouter(bool inside, Vec& first_douter, string& filename)
{
    int h = 0;
    ifstream ifs(filename.c_str(), ios::binary);
    if (!ifs || (h = GetColumn(ifs)) <= 0) return;
    if (!noout) cout << "#-Reading " << filename << endl;
    if (!noout) cout << "column size: " << h << endl;
    for (LEN i = 0; i <= _constraint; i++) {
        Vec temp = Vec();
        if (ReadBinVec(h, temp, ifs)) {
            if ((inside && i == 2*_constraint) || (!inside && i == 0)) {
                first_douter = temp;
                break;
            }
        } else  break;
    }
}

void ParasoR::ReadFirstDouter(bool inside, Vec& first_douter, string& filename)
{
    string str;
    ifstream ifs(filename.c_str());
    if (!ifs) return;
    if (!noout) cout << "#-Reading " << filename << endl;
    for (LEN i = 0; getline(ifs, str); i++) {
        Vec temp = Vec();
        ReadVec(temp, str);
        if ((inside && i == 2*_constraint) || (!inside && i == 0)) {
            first_douter = temp;
            break;
        }
    }
}

void ParasoR::ReadFirstDouter(bool inside)
{
    Vec first_douter;
    string filename = (inside) ? GetShrunkFileList(File::Shrunk, inside, id-1)
                                : GetShrunkFileList(File::Shrunk, inside, id+1);
    if ((inside && id > 0) || (!inside && id < chunk-1)) {
        (binary) ? ReadBinFirstDouter(inside, first_douter, filename) : ReadFirstDouter(inside, first_douter, filename);
        if (inside) alpha.douter[0] = first_douter;
        else beta.douter[static_cast<LEN>(beta.douter.size())-1] = first_douter;
    }
}

void ParasoR::ConnectDoSaved(bool inside, Vec& old_douter)
{
    if (inside) {
        ReadFirstDouter(true);
        ConnectInDo(old_douter, alpha.douter, id, GetShrunkFileList(File::Whole, true, id), false);
    } else {
        ReadFirstDouter(false);
        ConnectOutDo(old_douter, beta.douter, id, GetShrunkFileList(File::Whole, false, id), false);
    }
}

bool ParasoR::ConnectSavedFiles(bool keep_flag, bool inside)
{
    string outfile = GetDoFile(inside);
    bool flag = true;
    for (int i = 0; i < chunk; i++) {
        string file = GetShrunkFileList(File::Whole, inside, i);
        bool error = (binary) ? ReadBinStemToSingleFile(file, outfile, i != 0)
                 : ReadStemToSingleFile(file, outfile, i != 0, "");
        if (error) return false;
    }
    if (!binary)
        TabToNewLine(outfile);
    if (!noout) cout << "#-Complete " << GetDoFile(inside) << endl;
    if (!keep_flag) RemoveTempFilesSaveMem(inside);
    return true;
}

bool ParasoR::ConnectSavedFiles(bool keep_flag)
{
    bool flagi = ConnectSavedFiles(keep_flag, true);
    bool flago = ConnectSavedFiles(keep_flag, false);
    return flagi || flago;
}

/* ///////////////////////////////////////////// */




DOUBLE ParasoR::bpp(LEN i, LEN j, bool deb)
{
    if (i > j) swap(i, j);
    if (std::abs(i-j) < TURN) return 0.0;
    DOUBLE stack = Logsum(Stem(alpha, i, j-1), Stem(beta, i-1, j), LogLoopEnergy(i, j, i+1, j-1, seq));
    DOUBLE stemend = Logsum(Stemend(alpha, i, j-1), Stem(beta, i-1, j));
    DOUBLE temp = 0.0;
    if (!Is_INF(stack)) temp += exp(Logsum(stack, -Outer(alpha, seq.length)));
    if (!Is_INF(stemend)) temp += exp(Logsum(stemend, -Outer(alpha, seq.length)));
    if (deb) {
        cout << "stack " << Stem(alpha, i, j-1) << " " << Stem(beta, i-1, j) << " " << stack << endl;
        cout << "stemend " << Stemend(alpha, i, j-1) << " " << Stem(beta, i-1, j) << " " << stemend << endl;
        cout << "dpp " << i << " " << j << " " << temp << endl;
        cout << "Catch error: caused by the irregular outer file." << endl;
        exit(1);
    }
    return temp;
 }

DOUBLE ParasoR::bppDelta(LEN i, LEN j, bool deb)
{
    if (i > j) swap(i, j);
    if (std::abs(i-j) < TURN) return 0.0;
    DOUBLE stack = Logsum(Stem(alpha, i, j-1), Stem(beta, i-1, j), LogLoopEnergy(i, j, i+1, j-1, seq));
    DOUBLE stemend = Logsum(Stemend(alpha, i, j-1), Stem(beta, i-1, j));
    DOUBLE temp = exp(Logsumexp(stack, stemend));
    if (deb && temp-1.0 > 1e-6) {
        cout << "#----error " << i << " " << j << endl;
        cout << "#---- " << Stem(alpha, i, j-1) << " " << Stem(beta, i-1, j) << " " << LogLoopEnergy(i, j, i+1, j-1, seq) << endl;
        cout << "#---- " << Stemend(alpha, i, j-1) << " " << Stem(beta, i-1, j) << endl;
        cout << "dpp " << exp(Logsumexp(stack, stemend)) << endl;
        cout << "Catch error: caused by outer file." << endl;
        exit(1);
    }
    return temp;
}

void ParasoR::WriteBpp(Mat& data)
{
    data.clear();
    data = Mat(seq.length, Vec(_constraint, 0.0));
    for (LEN i = 1; i < seq.length; i++) {
        for (LEN j = i+1; j <= RightBpRange(i); j++)
            data[i-1][j-i-1] = bpp(i, j);
    }
}

void ParasoR::WriteAcc(Mat& data)
{
    data.clear();
    data = Mat(seq.length, Vec(_constraint+1, 0.0));
    for (LEN i = 1; i <= seq.length; i++) {
        for (LEN j = i; j <= min(seq.length, i+_constraint); j++) {
            data[i-1][j-i] = acc(i, j);
            if (data[i-1][j-i] > 1.0) {
                if (data[i-1][j-i]-1.0 >= 1.0e-11) {
                    ;//cerr << "error? " << i << " " << j << " " << data[i-1][j-i]-1.0 << endl;
                }
                data[i-1][j-i] = 1.0;
            }
        }
    }
}

void ParasoR::WriteStemProb(Vec& stem, const Mat& mat)
{
    stem.clear();
    stem = Vec(static_cast<LEN>(mat.size()));
    for (LEN i = 0; i < static_cast<LEN>(mat.size()); i++) {
        for (LEN j = 0; j < static_cast<LEN>(mat[i].size()); j++) {
            stem[i] += mat[i][j];
            stem[i+j+1] += mat[i][j];
        }
    }
}

void ParasoR::WriteStemProb(Vec& stem)
{
    stem.clear();
    stem = Vec(seq.length, 0.0);
    for (LEN i = 1; i <= seq.length; i++) {
        for (LEN j = i+1; j <= i+_constraint && j <= seq.length; j++) {
            DOUBLE value = bpp(i, j);
            stem[i-1] += value;
            stem[j-1] += value;
        }
    }
}

void ParasoR::OutputStemProb(bool _acc)
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    for (LEN i = 1; i <= seq.length; i++) {
        DOUBLE P[3] = { 1.0, 0.0, 0.0 };
        for (LEN j = max(static_cast<LEN>(1), i-_constraint); j <= min(seq.length, i+_constraint); j++) {
            if (!_acc && j == i) continue;
            DOUBLE value = ((_acc) ? acc(i, j) : bpp(i, j));
            if (_acc) { cout << "*\t" << i << "\t" << j << "\t" << value << endl; continue; }
            P[0] -= value;
            (j < i) ? P[2] += value : P[1] += value;
        }
        if (!_acc) cout << "*\t" << seq.strget(i) << "\t" << i << "\t" << P[1] << "\t" << P[2] << endl;
    }
}

void ParasoR::OutputBppCond(bool _acc)
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    if (_acc) {
        for (LEN i = 1; i <= seq.length; i++) {
            DOUBLE value = acc(i, i+1);
            cout << "*\t" << i << "\t" << value << endl;
        }
    } else {
        for (LEN i = 1; i <= seq.length; i++) {
            for (LEN j = max(static_cast<LEN>(1), i-_constraint); j <= min(seq.length, i+_constraint); j++) {
                if (j == i) continue;
                DOUBLE value = bpp(i, j);
                cout << "*\t" << i << "\t" << j << "\t" << value << endl;
            }
        }
    }
}

string ParasoR::GetDoFile(bool inside)
{
    ostringstream os;
    if (inside) {
        os << TMP << name << "douter_inside_" << _constraint;
        if (!binary) os << ".txt";
    } else {
        os << TMP << name << "douter_outside_" << _constraint;
        if (!binary) os << ".txt";
    }
    return os.str();
}

string ParasoR::GetStemFile(bool acc, bool prof, bool boundary)
{
    ostringstream os;
    if (acc && prof) {
        os << STMP << name << "prof_" << _window << "_" << _constraint;
    } else if (acc) {
        os << STMP << name << "acc_" << _window << "_" << _constraint;
    } else if (boundary) {
        os << STMP << name << "bound_" << _constraint;
    } else {
        os << STMP << name << "stem_" << _constraint;
    }
    if (!binary) os << ".txt";
    return os.str();
}

string ParasoR::GetDividedStemFile(bool acc, bool prof, bool boundary)
{
    ostringstream os;
    if (acc && prof) {
        os << STMP << name << "profd_" << id << "_" << chunk << "_" << _window << "_" << _constraint;
    } else if (acc) {
        os << STMP << name << "accd_" << id << "_" << chunk << "_" << _window << "_" << _constraint;
    } else if (boundary) {
        os << STMP << name << "bound_" << id << "_" << chunk << "_" << _constraint;
    } else {
        os << STMP << name << "stemd_" << id << "_" << chunk << "_" << _constraint;
    }
    if (!binary) os << ".txt";
    return os.str();
}

string ParasoR::GetIDFileList(string prefix, int tid)
{
    ostringstream os;
    os << TMP << name << prefix << tid << "_" << chunk << "_" << _constraint;
    if (!binary) os << ".txt";
    return os.str();
}

string ParasoR::GetShrunkFileList(int preID, bool inside, int tid)
{
    if (tid < 0) tid = id;
    string prefix;
    switch (preID) {
        case File::Part: prefix = "part"; tid = -1; break;
        case File::Shrunk: prefix = "shrunk"; break;
        case File::Whole: prefix = "whole"; break;
    }
    if (inside) {
        return GetIDFileList(prefix+"_inside_douter_", tid);
    } else {
        return GetIDFileList(prefix+"_outside_douter_", tid);
    }
}

string ParasoR::GetTempFileList(bool inside, int tid)
{
    if (tid < 0) tid = id;
    if (inside) {
        return GetIDFileList("temp_inside_douter_", tid);
    } else {
        return GetIDFileList("temp_outside_douter_", tid);
    }
}

void ParasoR::WriteDouterTemp(ofstream& ofs, Mat& mat)
{
    ostream_iterator<DOUBLE> oit(ofs, "\t");
    for (Mat::iterator it = mat.begin(); it != mat.end(); it++) {
        copy(it->begin(), it->begin()+(int)it->size(), oit);
        ofs << endl;
    }
}

void ParasoR::WriteBinDouterTemp(ofstream& ofs, Mat& mat)
{
    int h = mat[0].size();
    ofs.write((char*)&h, sizeof(int));
    for (LEN i = 0; i < static_cast<LEN>(mat.size()); i++) {
        for (LEN j = 0; j < static_cast<LEN>((mat[i].size())); j++) {
            ofs.write((char*)&(mat[i][j]), sizeof(DOUBLE));
        }
    }
}

void ParasoR::StoreDouterTemp(bool inside)
{
    string file = GetTempFileList(inside);
    if (binary) {
        ofstream ofs(file.c_str(), ios::binary);
        WriteBinDouterTemp(ofs, (inside) ? alpha.douter : beta.douter);
        ofs.close();
        if (!noout) cout << "#-Written (binary) " << file << endl;
    } else {
        ofstream ofs(file.c_str());
        ofs.precision(PREC);
        WriteDouterTemp(ofs, (inside) ? alpha.douter : beta.douter);
        if (!noout) cout << "#-Written " << file << endl;
        ofs.close();
    }
}

void ParasoR::LeaveDouterEdgeRegion(bool inside)
{
    if (inside) {
        Mat::const_iterator first = alpha.douter.begin()+max(static_cast<LEN>(0), static_cast<LEN>(alpha.douter.size())-2*_constraint-1);
        Mat::const_iterator last = alpha.douter.end();
        Mat temp(first, last);
        alpha.douter = temp;
    } else {
        Mat::const_iterator first = beta.douter.begin();
        Mat::const_iterator last = beta.douter.begin()+min(static_cast<LEN>(beta.douter.size()), 2*_constraint+1);
        Mat temp(first, last);
        beta.douter = temp;
    }
}

void ParasoR::StoreDouterTempShrunk(bool inside)
{
    string file = GetShrunkFileList(File::Shrunk, inside);
    bool flag = false;
    LeaveDouterEdgeRegion(inside);
    if (binary) {
        ofstream ofs(file.c_str(), ios::binary);
        WriteBinDouterTemp(ofs, (inside) ? alpha.douter : beta.douter);
        ofs.close();
        if (!noout) cout << "#-Written (binary) " << file << endl;
    } else {
        ofstream ofs(file.c_str());
        ofs.precision(PREC);
        WriteDouterTemp(ofs, (inside) ? alpha.douter : beta.douter);
        if (!noout) cout << "#-Written " << file << endl;
        ofs.close();
    }
}

void ParasoR::StoreDouter(string filename, Vec& douter, bool inside, bool app)
{
    if (douter.size() == 0) {
        cerr << "Cannot store douter (may have a problem in reading temp douter file)" << endl;
        abort();
    }
    LEN start = (inside) ? 1 : 0;
    LEN end = (inside) ? 0 : -1;
    if (binary) {
        int h = douter.size()+static_cast<int>(end-start);
        ofstream ofs(filename.c_str(), ((app) ? ios::app : ios::trunc) | ios::binary);
        ofs.write((char*)&h, sizeof(int));
        for (LEN i = start; i < static_cast<LEN>(douter.size())+end; i++)
            ofs.write((char*)(&douter[i]), sizeof(DOUBLE));
        ofs.close();
    } else {
        ofstream ofs(filename.c_str(), (app) ? ios::app : ios::trunc);
        ofs.precision(PREC);
        ostream_iterator<DOUBLE> oit(ofs, "\t");
        copy(douter.begin()+start, douter.begin()+(int)douter.size()+end, oit);
        ofs << endl;
        ofs.close();
    }
    if (!noout) cout << "#-Written " << filename << endl;
}

void ParasoR::StoreStem(string filename, Vec& vec, bool acc)
{
    assert(vec.size() > 0);
    if (acc) ;
    else if (*max_element(vec.begin(), vec.end())-1.0 >= 0.1 || *min_element(vec.begin(),vec.end()) <= -0.1)
        throw "Value overflow or underflow";
    if (binary) {
        LEN const_size = vec.size();
        ofstream ofs(filename.c_str(), ios::binary);
        ofs.write((char*)&const_size, sizeof(int));
        for (LEN i = 0; i < _end-_start && i < static_cast<LEN>(vec.size()); i++) {
            if (!acc) vec[i] = min(max((DOUBLE)0, vec[i]), (DOUBLE)1);
            ofs.write((char*)(&vec[i]), sizeof(DOUBLE));
        }
        ofs.close();
    } else {
        ofstream ofs(filename.c_str(), ios::trunc);
        ofs.precision(PREC);
        ostream_iterator<DOUBLE> oit(ofs, "\t");
        copy(vec.begin(), vec.begin()+(int)vec.size(), oit);
        ofs.close();
    }
}

// void ParasoR::StoreProf(string filename, vector<char>& vec)
// {
//     assert(vec.size() > 0);
//     if (binary) {
//         int const_size = vec.size();
//         ofstream ofs(filename.c_str(), ios::binary);
//         ofs.write((char*)&const_size, sizeof(int));
//         for (int i = 0; i < _end-_start && i < (int)vec.size(); i++) {
//             ofs.write((char*)(&vec[i]), sizeof(char));
//         }
//         ofs.close();
//     } else {
//         ofstream ofs(filename.c_str(), ios::trunc);
//         ostream_iterator<char> oit(ofs, "\t");
//         copy(vec.begin(), vec.begin()+(int)vec.size(), oit);
//         ofs.close();
//     }
// }

void ParasoR::PrintMat(bool is_alpha)
{
    if (is_alpha) {
        cout << "#--alpha" << endl;
        alpha.Print(seq.substr(alpha.istart));
    } else {
        cout << "#--beta" << endl;
        beta.Print(seq.substr(alpha.istart));
    }
}

void ParasoR::InitMatrix(LEN row, LEN col)
{
    alpha = Matrix(row, col, true);
    beta = Matrix(row, col, false);
}

void ParasoR::Init(bool full, bool connect)
{
    if (full) InitMatrix(_end-_start+1, _constraint);
    else {
        if (cut && !connect) {
            seq.CutSequence(max(static_cast<LEN>(0), _start-CONST*_constraint), min(seq.length, _end+CONST*_constraint));
        }
        InitMatrix(CONST*_constraint, _constraint);
    }
    SetIndex(_start, _end, true);
}

void ParasoR::InitBpp(bool full, bool set)
{
    if (!alpha.isSet() || !beta.isSet()) {
        if (full) InitMatrix(_end-_start+1, _constraint);
        else InitMatrix(CONST*_constraint, _constraint);
    }
    LEN left = max(_start-CONST*_constraint, static_cast<LEN>(0));
    LEN right = min(_end+CONST*_constraint, seq.length);
    if (!full && cut && set) {
        seq.CutSequence(left, right);
    }
    SetIndex(left, right, false, set);
}

void ParasoR::SetSequence(const string& filename, int seqID, LEN length)
{
    seq = Sequence(filename, seqID, length);
}

void ParasoR::SetSequence(const string& sequence, bool out)
{
    if (out) {
        cout << "# " << sequence.substr(0, 50);
        if (sequence.length() > 50) cout << "...";
        cout << endl;
    }
    string str = "$" + sequence;
    seq = Sequence(str);
}

void ParasoR::RemoveTempFilesSaveMem(bool inside)
{
    string file;
    for (int prefix = 0; prefix < 3; prefix++) {
        for (int i = 0; i < chunk; i++) {
            if (memory) {
                file = GetShrunkFileList(prefix, inside, i);
            } else
                file = GetTempFileList(inside, i);
            remove(file.c_str());
            if (!noout && i == 0)
                cout << "#-Remove " << file << " and others..." << endl;
        }
    }
}

void ParasoR::RemoveTempFiles(bool inside, int prefix)
{
    cout << "#-Remove Temp files." << endl;
    string file;
    for (int i = 0; i < chunk; i++) {
        if (memory) {
            file = GetShrunkFileList(prefix, inside, i);
        } else
            file = GetTempFileList(inside, i);
        remove(file.c_str());
    }
}

void ParasoR::RemoveStem(bool acc, bool prof)
{
    cout << "#-Remove Stem files" << endl;
    for (int i = 0; i < chunk; i++) {
        ostringstream os;
        id = i;
        string file = GetDividedStemFile(acc, prof);
        remove(file.c_str());
    }
}

void ParasoR::ConcatDo()
{
    if (binary) {
        Douter_concat::ConcatBin(GetDoFile(true));
        Douter_concat::ConcatBinReverse(GetDoFile(false));
        if (!noout) cout << "#-Complete Do file." << endl;
    } else {
        Douter_concat::Concat(GetDoFile(true));
        Douter_concat::ConcatReverse(GetDoFile(false));
        if (!noout) cout << "#-Complete Do file." << endl;
    }
}

void ParasoR::DivideChunk(Arg& arg, bool shrink)
{
    if (arg.str.length() == 0 && arg.length == 0) return;
    ParasoR rfold;
    rfold.SetBasicParam(arg, shrink);
    if (rfold.seq.length/arg.chunk > static_cast<LEN>(INT_MAX)) {
        cout << INT_MAX << " " << rfold.seq.length << endl;
        cerr << "Too long sequence for a single chunk.\nPlease increase the chunk number more than  " << rfold.seq.length/static_cast<LEN>(INT_MAX) << endl;
    } else if (rfold.SetChunkId(arg.id, arg.chunk)) {
        if (rfold.seq.length/rfold.chunk) {
            rfold.cut = true;
            rfold.CalcDeltaInOut();
        } else if (rfold.id == 0 && rfold.chunk == 1) {
            Connect(rfold, shrink);
        }
    } else {
        cerr << "Too many chunk for this sequence len: " << rfold.seq.length << endl;
    }
}

void ParasoR::PreviousCalculation(Arg& arg, bool shrink)
{
    if (arg.str.length() == 0) return;
    ParasoR rfold;
    rfold.SetWindow(arg.window, arg.acc_flag & !arg.prof_flag);
    rfold.SetBasicParam(arg, shrink);
    rfold.SetRange(1, arg.str.length());
    if (arg.prof_flag) {
        if (arg.acc_flag) {
            if (arg.mea_flag) rfold.CalcAllAtOnce(Out::PROFIM, arg.image, arg.gamma);
            else rfold.CalcAllAtOnce(Out::PROF);
        } else
            rfold.CalcAllAtOnce(Out::MOTIF);
    } else if (arg.acc_flag) {
        rfold.CalcAllAtOnce(Out::ACC);
    } else if (arg.mea_flag) {
        rfold.CalcAllAtOnce(Out::BPPIM, arg.image, arg.gamma);
    } else if (arg.stem_flag) {
        rfold.CalcAllAtOnce(Out::STEM);
    } else if (arg.entro_flag) {
        rfold.CalcAllAtOnce(Out::ENTRO);
    } else {
        rfold.CalcAllAtOnce(Out::BPP, arg.image, arg.minp);
    }
}

void ParasoR::main(Arg& arg, bool shrink)
{
    if (arg.str.length() == 0 && arg.length == 0) return;
    ParasoR rfold;
    rfold.SetWindow(arg.window, arg.acc_flag & !arg.prof_flag);
    if (arg.mtype == Arg::Mut::Ins && arg.constraint == arg.str.length())
        rfold.SetBasicParam(arg, false);
    else
        rfold.SetBasicParam(arg, shrink);
        if (arg.end < 0) {
            arg.end = max(arg.length, static_cast<LEN>(arg.str.length()));
            if (arg.start < 0) arg.start = 0;
        }
    if (!noout) cout << "#-Calculate from " << arg.start << " to " << arg.end << endl;
    if (arg.end > arg.start) {
        rfold.SetBpRange(arg.start, arg.end);
        rfold.cut = true;
        if (arg.mtype >= 0) {
            rfold.MutatedStem(arg);
        } else if (arg.mea_flag) {
            if (arg.prof_flag)
                rfold.CalcProf(arg.acc_flag);
            rfold.CalcMEA(true, arg.image, arg.prof_flag);
        } else if (arg.prof_flag) {
            rfold.CalcProf(arg.acc_flag);
        } else if (arg.acc_flag) {
            rfold.SetWindow(max(2, arg.window));
            rfold.CalcAcc();
        } else if (arg.stem_flag) {
            rfold.CalcStem(false, arg.boundary);
        } else if (arg.entro_flag) {
            rfold.CalcEntropy(true);
        } else {
            rfold.CalcBpp(true, arg.minp);
        }
    }

}


}
