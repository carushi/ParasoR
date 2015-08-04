#include "part_func.hh"

namespace Rfold {

DOUBLE ParasoR::HairpinAcc(LEN x1, LEN x2)
{
    DOUBLE value = -INF;
    for (LEN i = LeftRange(x2); i+1 < x1; i++)
        for (LEN j = max(x2+1, i+TURN+2); j <= RightRange(i); j++) {
            if (IsPair(seq.slidebp(i+1, j)))
                value = Logsumexp(value, Logsum(Stem(beta, i, j), LogHairpinEnergy(i+1, j, seq)));
        }
    return value;
}

DOUBLE ParasoR::InnerLoopAcc(LEN i, LEN j, LEN k, LEN l, bool deb) {
    if (deb) {
        cout << k+1-i << "-----" << j+1-l << endl;
        cout << i << " " << j << " " << k << " " << l << endl;
    }
    return Logsum(Stemend(beta, i, j), LogLoopEnergy(i, j+1, k+1, l, seq), Stem(alpha, k, l));
}

DOUBLE ParasoR::InteriorAcc(LEN x1, LEN x2, bool exbulge, bool exinter)
{
    DOUBLE value = -INF;
    for (LEN i = LeftBpRange(x2+TURN+1); i < x1; i++) { /* i -- x1 -- x2 -- k loop l -- j */
        for (LEN j = x2+TURN+1; j+1 <= RightBpRange(i); j++) {
            if (!IsPair(seq.slidebp(i, j+1))) continue;
            for (LEN k = x2; k+TURN+1 < j && k-i <= MAXLOOP; k++) {
                if (!exbulge) value = Logsumexp(value, InnerLoopAcc(i, j, k, j)); /* Bulge */
                if (!exinter) {
                    for (LEN l = max(k+TURN+2, j-MAXLOOP+(k-i)); l <= j-1; l++) {
                        value = Logsumexp(value, InnerLoopAcc(i, j, k, l));
                    }
                }
            }
        }
    }
    for (LEN i = LeftBpRange(x2); i < x1-TURN-1; i++) { /* i -- k loop l x1 -- x2 j */
         for (LEN j = x2; j+1 <= RightBpRange(i); j++) {
            if (!IsPair(seq.slidebp(i, j+1))) continue;
            if (!exbulge) {
                for (LEN l = max(i+TURN+2, j-MAXLOOP); l+1 <= x1; l++)
                    value = Logsumexp(value, InnerLoopAcc(i, j, i, l)); /* Bulge */
            }
            if (!exinter) {
                for (LEN k = i+1; k+TURN+2 < x1; k++) {
                    for (LEN l = max(k+TURN+2, j-MAXLOOP+(k-i)); l+1 <= x1; l++) {
                        if (IsPair(seq.slidebp(k+1, l)))
                            value = Logsumexp(value, InnerLoopAcc(i, j, k, l));
                    }
                }
            }
        }
    }
    return value;
}

DOUBLE ParasoR::BulgeAcc(LEN x1, LEN x2) {
    return InteriorAcc(x1, x2, false, true);
}

DOUBLE ParasoR::MultiAcc(LEN x1, LEN x2)
{
    DOUBLE value = -INF;
    for (LEN j = x2+1; j <= RightBpRange(x1); j++)
        value = Logsumexp(value, Logsum(Multi(beta, x1-1, j), (x2-x1+1)*logML_BASE, Multi(alpha, x2, j)));
    return value;
}

DOUBLE ParasoR::Multi2Acc(LEN x1, LEN x2)
{
    DOUBLE value = -INF;
    for (LEN i = LeftRange(x2); i < x1-1; i++)
        value = Logsumexp(value, Logsum(Multi2(beta, i, x2), (x2-x1+1)*logML_BASE, Multi2(alpha, i, x1-1)));
    return value;
}

DOUBLE ParasoR::acc(LEN i, LEN j)
{
    if (i > j) swap(i, j);
    DOUBLE value = HairpinAcc(i, j);
    value = Logsumexp(value, InteriorAcc(i, j, false));
    value = Logsumexp(value, MultiAcc(i, j));
    value = Logsumexp(value, Multi2Acc(i, j));
    value = Logsumexp(value, Logsum(beta.outer[j], alpha.outer[i-1]));
    value = Logsum(value, -alpha.outer[seq.length]);
    return ((ene) ? Energy(value) : exp(value));
}

DOUBLE ParasoR::accDelta(LEN i, LEN j, DOUBLE uij)
{
    if (i > j) swap(i, j);
    DOUBLE value = HairpinAcc(i, j);
    value = Logsumexp(value, InteriorAcc(i, j, false));
    value = Logsumexp(value, MultiAcc(i, j));
    value = Logsumexp(value, Multi2Acc(i, j));
    value = Logsumexp(value, uij);
    return ((ene) ? Energy(value) : exp(value));
}

DOUBLE ParasoR::Hairpinprof(LEN x) {
    return HairpinAcc(x, x);
}

DOUBLE ParasoR::Stemprof(LEN i) {
    DOUBLE value = 0.;
    for (LEN j = LeftBpRange(i); j <= i-TURN; j++)
        value +=  (delta) ? bppDelta(j, i) : bpp(j, i);
    for (LEN j = RightBpRange(i); j >= i+TURN; j--)
        value +=  (delta) ? bppDelta(i, j) : bpp(i, j);
    return value;
}

DOUBLE ParasoR::Multiprof(LEN i) {
    return Logsumexp(MultiAcc(i, i), Multi2Acc(i, i));
}

DOUBLE ParasoR::Interiorprof(LEN i) {
    return InteriorAcc(i, i, true);
}

DOUBLE ParasoR::Bulgeprof(LEN i) {
    return BulgeAcc(i, i);
}

DOUBLE ParasoR::Otherprof(Vec& profile) {
    return 1.-accumulate(profile.end()-TYPE, profile.end(), 0.);
}

void ParasoR::profile(LEN i, Vec& profvec)
{
    DOUBLE Z = alpha.outer[seq.length];
    profvec.push_back(exp(Logsum(Bulgeprof(i), -Z)));
    profvec.push_back(exp(alpha.outer[i-1]+beta.outer[i]-alpha.outer[seq.length]));
    profvec.push_back(exp(Logsum(Hairpinprof(i), -Z)));
    profvec.push_back(exp(Logsum(Multiprof(i), -Z)));
    profvec.push_back(Stemprof(i));
    profvec.push_back(exp(Logsum(Interiorprof(i), -Z)));
    if (ddebug) {
        profvec.push_back(Otherprof(profvec));
        profvec.push_back(1.-(Stemprof(i)+acc(i, i)));
    }
}

void ParasoR::profileDelta(LEN i, DOUBLE uij, Vec& profvec)
{
    profvec.push_back(exp(Bulgeprof(i)));
    profvec.push_back(exp(uij));
    profvec.push_back(exp(Hairpinprof(i)));
    profvec.push_back(exp(Multiprof(i)));
    profvec.push_back(Stemprof(i));
    profvec.push_back(exp(Interiorprof(i)));
    if (ddebug) {
        profvec.push_back(Otherprof(profvec));
    }
}

void ParasoR::profileDeltaLinear(LEN i, LEN j, DOUBLE uij, Vec& profvec)
{
    profvec = Vec(TYPE, 0.);
    profvec[4] = StemLinear(i, j);
    if (i == j)
        profvec[1] = expInf(uij);
    if (i == j)
        profvec[3] = expInf(MultiLinear(i));
    if (i <= 1 || j >= seq.length) return;
    profvec[0] = expInf(BulgeLinear(i, j));
    profvec[2] = expInf(HairpinLinear(i, j));
    profvec[5] = expInf(InteriorLinear(i, j));
}

DOUBLE ParasoR::StemLinear(LEN i, LEN j) {
    return (delta) ? bppDelta(i, j) : bpp(i, j);
}

DOUBLE ParasoR::InteriorLinearLeft(LEN x1, LEN x2, LEN i, LEN k, bool exbulge, bool exinter)
{
    DOUBLE value = -INF;
    if (i < LeftBpRange(x2+TURN+1)) return value; /* i -- x1 -- x2 -- k loop l -- j */
    for (LEN j = x2+TURN+1; j+1 <= RightBpRange(i); j++) {
        if (!IsPair(seq.slidebp(i, j+1))) continue;
        if (k+TURN+1 >= j || k-i > MAXLOOP) continue;
        if (!exbulge) value = Logsumexp(value, InnerLoopAcc(i, j, k, j)); /* Bulge */
        if (!exinter) {
            for (LEN l = max(k+TURN+2, j-MAXLOOP+(k-i)); l <= j-1; l++) {
                value = Logsumexp(value, InnerLoopAcc(i, j, k, l));
            }
        }
    }
    return value;
}

DOUBLE ParasoR::InteriorLinearRight(LEN x1, LEN x2, LEN j, LEN l, bool exbulge, bool exinter)
{
    DOUBLE value = -INF;
    for (LEN i = LeftBpRange(x2); i < x1-TURN-1; i++) { /* i -- k loop l x1 -- x2 j */
        if (j+1 > RightBpRange(i)) continue;
        if (!IsPair(seq.slidebp(i, j+1))) continue;
        if (!exbulge) {
            if (l >= max(i+TURN+2, j-MAXLOOP))
                value = Logsumexp(value, InnerLoopAcc(i, j, i, l)); /* Bulge */
        }
        if (!exinter) {
            for (LEN k = i+1; k+TURN+2 < x1; k++) {
                if (l < max(k+TURN+2, j-MAXLOOP+(k-i))) continue;
                if (IsPair(seq.slidebp(k+1, l)))
                    value = Logsumexp(value, InnerLoopAcc(i, j, k, l));
            }
        }
    }
    return value;
}

DOUBLE ParasoR::InteriorLinear(LEN x1, LEN x2, bool exbulge, bool exinter)
{
    return Logsumexp(InteriorLinearLeft(x1, x2, x1-1, x2, exbulge, exinter),
                     InteriorLinearRight(x1, x2, x2, x1-1, exbulge, exinter));
}

DOUBLE ParasoR::BulgeLinear(LEN x1, LEN x2) {
    return InteriorLinear(x1, x2, false, true);
}

DOUBLE ParasoR::InteriorLinear(LEN x1, LEN x2) {
    return InteriorLinear(x1, x2, true, false);
}

DOUBLE ParasoR::HairpinLinear(LEN x1, LEN x2)
{
    LEN i = x1-2, j = x2+1;
    if (j < i+TURN+2 || j > RightRange(i)) return -INF;
    if (!IsPair(seq.slidebp(i+1, j))) return -INF;
    return Logsum(Stem(beta, i, j), LogHairpinEnergy(i+1, j, seq));
}

DOUBLE ParasoR::MultiLinear(LEN x)
{
    DOUBLE value = -INF;
    for (LEN j = x+1; j <= RightBpRange(x); j++)
        value = Logsumexp(value, Logsum(Multi(beta, x-1, j), logML_BASE, Multi(alpha, x, j)));
    for (LEN i = LeftRange(x); i < x-1; i++)
        value = Logsumexp(value, Logsum(Multi2(beta, i, x), logML_BASE, Multi2(alpha, i, x-1)));
    return value;
}

void ParasoR::AddCumProf(Vec& cump, LEN pos, Vec& prof)
{
    LEN index = pos*TYPE;
    transform(prof.begin()+index, prof.begin()+index+TYPE, cump.begin(), prof.begin()+index, plus<DOUBLE>());
}

void ParasoR::GetProfsLinearRight(LEN pos, Vec& prof, LEN bstart, LEN end)
{
    for (LEN i = LeftBpRange(pos); i <= end; i++) {
        Vec temp;
        profileDeltaLinear(i, pos, -INF, temp);
        for (LEN h = i; h <= end; h++) {
            if (h >= bstart && h <= end) AddCumProf(temp, h-bstart, prof);
            temp[4] = 0.; // Stem profile is not regional.
        }
    }
}

void ParasoR::GetProfsLinear(LEN pos, Vec& prof, LEN bstart, LEN end, DOUBLE uxx)
{
    DOUBLE uij = ExpandLocalOuter(uxx, pos, pos-1, pos);
    for (LEN i = LeftBpRange(pos); i <= pos; i++) {
        Vec temp;
        profileDeltaLinear(i, pos, uij, temp);
        if (pos <= end)
            AddCumProf(temp, pos-bstart, prof);
        temp[1] = temp[3] = 0.;
        for (LEN h = i; h < pos; h++) {
            if (h >= bstart && h <= end) AddCumProf(temp, h-bstart, prof);
            temp[4] = 0.; // Stem profile is not regional.
        }
    }
}

void ParasoR::GetProfs(LEN pos, Vec& prof, LEN bstart, DOUBLE uij, bool store)
{
    if (delta) {
        Vec temp;
        profileDelta(pos, uij, temp);
        AddCumProf(temp, pos-bstart, prof);
    }
    else {
        if (!store) prof = Vec();
        profile(pos, prof);
    }
}

}