#include "part_func.hh"
/**
 * Hard constraint and MFE structure.
 **/

namespace Rfold {

bool ParasoR::SameSequenceHC(LEN i, LEN j)
{
    if (i > j) return false;
    if (!hard) return true;
    if (seq.hcsubstr(i, j-i+1).find('p') != string::npos) return false;
    // if (seq.hcsubstr(i, j-i+1).find('$') != string::npos) {
        // for (LEN k = i; k <= j && hard; k++) {
        //     if (seq.hcget(k) != seq.hcget(i)) {
        //         cout << "break" << endl;
        //         return false;
        //     }
        // }
    // } else {
    //     for (LEN k = i; k <= j && hard; k++) {
    //         if (!IsLoopHC(k)) {
    //             cout << "break" << endl;
    //             return false;
    //         }
    //     }
    // }
    return true;
}

bool ParasoR::InnerLoopHC(LEN i, LEN j, LEN p, LEN q) {
    if (!hard) return true;
    else if (p < i || q > j) return false;
    else if (seq.hcsubstr(i, p-i+1).find('p') != string::npos || seq.hcsubstr(i, p-i+1).find('$') != string::npos) {
        return false;
    } else if (seq.hcsubstr(q, j-q+1).find('p') != string::npos || seq.hcsubstr(q, j-q+1).find('$') != string::npos) {
        return false;
    }
    return true;
}

bool ParasoR::IsLoopHC(LEN i) {
    return !hard || (seq.hcget(i) != 'p' && seq.hcget(i) != '$');
}

bool ParasoR::IsPairHC(LEN i, LEN j) {
    return IsPair(seq.slidebp(i, j)) && (!hard || (seq.hcget(i) != 'l' && seq.hcget(j) != 'l' && seq.hcget(i) != '$' && seq.hcget(j) != '$'));
}

DOUBLE ParasoR::GetInStemHC(LEN i, LEN j)
{
    int bp1 = seq.slidebp(i+1, j);
    DOUBLE temp = -INF;
    if (IsPairHC(i+1, j)) {
        temp = min_or_sum(temp, Stemend(alpha, i+1, j-1));
        if (!Is_INF(Stem(alpha, i+1, j-1)) && IsPairHC(i+2, j-1)) {
            temp = min_or_sum(temp, Logsum(Stem(alpha, i+1, j-1), LogLoopEnergy(i+1, j, i+2, j-1, seq)));
        }
    }
    return temp;
}

DOUBLE ParasoR::GetInMultiBifHC(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    for (LEN k = i+TURN+1; k < j-TURN; k++) {
        temp = min_or_sum(temp, Logsum(Multi1(alpha, i, k), Multi2(alpha, k, j)));
    }
    return temp;
}

DOUBLE ParasoR::GetInMulti2HC(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    if (IsLoopHC(j))
        temp = Logsum(Multi2(alpha, i, j-1), logML_BASE);
    if (!Is_INF(Stem(alpha, i, j))) {
        temp = min_or_sum(temp, Logsum(Stem(alpha, i, j), SumExtML(seq.slidebp(i+1, j), i, j+1, false, seq), logMLintern));
    }
    return temp;
}

DOUBLE ParasoR::GetInMulti1HC(LEN i, LEN j) {
    DOUBLE temp = min_or_sum(Multi2(alpha, i, j), Multibif(alpha, i, j));
    return temp;
}

DOUBLE ParasoR::GetInMultiHC(LEN i, LEN j) {
    DOUBLE temp = Multibif(alpha, i, j);
    if (IsLoopHC(i+1))
        temp = min_or_sum(temp, Logsum(Multi(alpha, i+1, j), logML_BASE));
    return temp;
}


DOUBLE ParasoR::GetInStemendHC(LEN i, LEN j)
{
    LEN p = i, q = j+1;
    DOUBLE temp = -INF;
    if (q > p+_constraint || q > seq.length || p <= 0 || !IsPairHC(p, q)) return temp;
    if (SameSequenceHC(p+1, q-1)) {
        if (hard && seq.hcsubstr(p+1, q-p-1).find('$') != string::npos) {// For interaction energy computation;
            if (seq.hcget(p+1) == '$' && seq.hcget(q-1) == '$')
                temp = 0.;
            else
                temp = SumExtML(seq.sliderbp(p, q), j, i+1, false, seq);
        } else {
            temp = LogHairpinEnergy(p, q, seq);
        }
    }
    for (LEN ip = i; ip < j-TURN-1; ip++) {
        for (LEN jp = max(j-MAXLOOP+(ip-i), ip+TURN+2); jp <= j && (j-jp)+(ip-i) > 0; jp++) {
            if (!Is_INF(Stem(alpha, ip, jp)) && InnerLoopHC(p+1, q-1, ip, jp+1)) {
                temp = min_or_sum(temp, Logsum(Stem(alpha, ip, jp), LogLoopEnergy(p, q, ip+1, jp, seq)));
            }
        }
    }
    temp = min_or_sum(temp, Logsum(Multi(alpha, i, j), SumExtML(seq.sliderbp(p, q), j, i+1, false, seq),
                                   logMLclosing+logMLintern));
    return temp;
}

/* ///////////////////////////////////////////// */

DOUBLE ParasoR::GetOutStemendHC(LEN i, LEN j) {
    if (!IsOnlyRange(i-1, j+1)) return -INF;
    else
        return Stem(beta, i-1, j+1);
}

DOUBLE ParasoR::GetOutMultiHC(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    if (IsOnlyRange(i-1, j)) {
        if (IsLoopHC(i)) {
            temp = Logsum(Multi(beta, i-1, j), logML_BASE);
        }
    }
    if (IsOnlyRange(i, j+1))
        temp = min_or_sum(temp, Logsum(Stemend(beta, i, j),
                                      SumExtML(seq.sliderbp(i, j+1), j, i+1, false, seq),
                                      logMLclosing+logMLintern));
    return temp;
}

DOUBLE ParasoR::GetOutMulti1HC(LEN i, LEN j)
{
    DOUBLE temp = -INF;
    for (LEN k = j+TURN+2; k < min(seq.length, i+_constraint+2); k++)
        temp = min_or_sum(temp, Logsum(Multibif(beta, i, k), Multi2(alpha, j, k)));
    return temp;
}

DOUBLE ParasoR::GetOutMulti2HC(LEN i, LEN j)
{
    DOUBLE temp = Multi1(beta, i, j);
    if (IsOnlyRange(i, j+1)) {
        if (IsLoopHC(j+1)) {
            temp = min_or_sum(temp, Logsum(Multi2(beta, i, j+1), logML_BASE));
        }
    }
    for (LEN k = max(static_cast<LEN>(0), j-_constraint-1); k < i-TURN-1; k++)
        temp = min_or_sum(temp, Logsum(Multibif(beta, k, j), Multi1(alpha, k, i)));
    return temp;
}

DOUBLE ParasoR::GetOutMultiBifHC(LEN i, LEN j) {
    return min_or_sum(Multi1(beta, i, j), Multi(beta, i, j));
}

DOUBLE ParasoR::GetOutStemHC(LEN i, LEN j) {
    return GetOutStemHC(i, j, Logsum(Outer(alpha, i), Outer(beta, j)));
}

DOUBLE ParasoR::GetOutStemHC(LEN i, LEN j, DOUBLE value)
{
    LEN p = i+1, q = j, type = seq.slidebp(p, q);
    DOUBLE temp = -INF;
    if (IsPairHC(p, q)) {
        temp = Logsum(value, SumExtML(type, i, j+1, true, seq));
        if (q-p-1 >= TURN) {
            for (LEN ip = i; ip >= LeftRange(j); ip--) {
                for (LEN jp = min(ip+_constraint-1, min(j+MAXLOOP-(i-ip), seq.length-1)); jp >= j && (i-ip)+(jp-j) > 0; jp--) {
                    if (IsPairHC(ip, jp+1) && InnerLoopHC(ip+1, jp, p-1, q+1)) {
                        temp = min_or_sum(temp, Logsum(Stemend(beta, ip, jp), LogLoopEnergy(ip, jp+1, p, q, seq)));
                    }
                }
            }
        }
        temp = min_or_sum(temp, Logsum(Multi2(beta, i, j), SumExtML(type, i, j+1, false, seq), logMLintern));
        if (IsOnlyRange(i-1, j+1))
            temp = min_or_sum(temp, Logsum(Stem(beta, i-1, j+1), LogLoopEnergy(i, j+1, p, q, seq)));
    }
    return temp;
}

void ParasoR::CalcInsideOuterHC(LEN j)
{
    assert(delta == false);
    if (j == 0) Outer(alpha, 0) = 0.0;
    else {
        DOUBLE temp = -INF;
        for (LEN k = LeftRange(j); k < j; k++) {
            // cout << j << " " << k << " " << temp << endl;
            temp = min_or_sum(temp, Logsum(Outer(alpha, k), Stem(alpha, k, j), SumExtML(seq.slidebp(k+1, j), k, j+1, true, seq)));
        }
        if (IsLoopHC(j)) {
            // cout << j << "loop" << endl;
            temp = min_or_sum(temp, Outer(alpha, j-1));
        }
        Outer(alpha, j) = temp;
    }
}

void ParasoR::CalcOutsideOuterHC(LEN j)
{
    assert(delta == false);
    if (j == seq.length) Outer(beta, seq.length) = 0.0;
    else {
        DOUBLE temp = -INF;
        for (LEN k = j+TURN; k <= RightRange(j); k++) {
            temp = min_or_sum(temp, Logsum(Stem(alpha, j, k), Outer(beta, k), SumExtML(seq.slidebp(j+1, k), j, k+1, true, seq)));
        }
        if (IsLoopHC(j+1)) {
            temp = min_or_sum(temp, Outer(beta, j+1));
        }
        Outer(beta, j) = temp;
    }
}

void ParasoR::CalcOuterHC()
{
    for (LEN j = 0; j <= seq.length; j++) {
        if (j >= TURN-1)
            CalcInside(j);
        CalcInsideOuterHC(j);
    }
    for (LEN i = seq.length; i >= 0; i--) {
        if (i <= seq.length-TURN)
            CalcInsideFromRight(i);
        CalcOutsideOuterHC(i);
    }
}


}