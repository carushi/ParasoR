#ifndef _SEQUENCE_HH
#define _SEQUENCE_HH

#include "pair_mat.hh"
#include "energy_const.hh"
#include "matrix.hh"
#include <cctype>

#include <vector>
#include <string>

namespace Rfold {
namespace Parameter {

using std::vector;
using std::string;
using std::min;
using std::max;

static char ComplementRNA(char c)
{
    if (c == '$') return c;
    c = toupper(c);
    switch(c) {
        case 'T': case 'U': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        case 'A': return 'U';
        case 'N': return 'N';
        default:
            return 'N';
        /*
        case 'W': case 'K': case 'Y': case 'H': case 'B': case 'D': case 'N': return 'A';
        case 'R': case 'M': case 'V': return 'U';
        case 'S': return 'C';
        */
    }
}

// static char CapitalRNA(char c)
// {
//     if (c == '$') return c;
//     c = toupper(c);
//     switch(c) {
//         case 'A': case 'C': case 'G': case'U': case 'N': return c;
//         case 'T': case 'K': case 'Y': case 'B': return 'U';
//         case 'W': case 'R': case 'M':
//         case 'H': case 'V': case 'D': return 'A';
//         case 'S': return 'C';
//         default:
//             std::cerr << "ambiguous character " << c << endl;
//             return 'N';
//     }
// }

static char CapitalRNA(char c)
{
    if (c == '$') return c;
    c = toupper(c);
    switch(c) {
        case 'A': case 'C': case 'G': case'U': case 'N': return c;
        case 'T': return 'U';
        default:
            return 'N';
    }
}

static void GetCapitalRNA(string& seq)
{
    string str;
    transform(seq.begin(), seq.end(), back_inserter(str), CapitalRNA);
    seq = str;
}

static void GetCompCapitalRNA(string& seq)
{
    string str;
    transform(seq.begin(), seq.end(), back_inserter(str), ComplementRNA);
    reverse(str.begin(), str.end());
    seq = str;
}
static void EraseWindowsNewline(string& seq)
{
    if (seq[seq.length()-1] == '\r')
        seq.erase(seq.begin()+seq.length()-1);
}

/**
 * @return Whether 'type' of base pair is possible.
 */
inline static bool IsPair(int type) {
    return (type > 0);
}
/**
 * @param sequence vector of base index.
 * @return type of base pair between 'i' and 'j'.
 */
inline static int bp(LEN i, LEN j, const vector<char>& sequence) {
    return BP_pair[(int)sequence[i]][(int)sequence[j]];
}
/**
 * @param Sequence vector of base index.
 * @return Reversed type of base pair between 'i' and 'j'.
 */
inline static int rbp(LEN i, LEN j, const vector<char>& sequence) {
    return rtype[bp(i, j, sequence)];
}

/**
 * A Sequence class.
 * Stores (fragmented) base sequence.
 */
class Sequence {
private:
    bool part;  // Don't include full sequence;
    bool lazy;  // Lazy evaluation for sequence import;
    LEN _shift; // index to be slided;
    string str; // base string;
    string filename;
    int seqID;
    vector<char> sequence;

    void CutVector(LEN start, LEN end) {
        sequence = vector<char>(sequence.begin()+start, sequence.begin()+end);
    }
    void CutString(LEN start, LEN end) {
        cout << "#-Cut " << start << " " << end << " " << str.length() << endl;
        str = str.substr(start, end-start);
    }
    /**
     * Set partial str.
     */
    void ReadPartialSeq(LEN start, LEN end)
    {
        int count = 0;
        LEN readthrough = 1;
        string temp;
        ifstream ifs(filename.c_str());
        str = (start == 0) ? "$" : "";
        while (getline(ifs, temp)) {
            if (temp.length() == 0) continue;
            EraseWindowsNewline(temp);
            if (temp[0] == '>') {
                count++;
                if (count > seqID) break;
            } else if (count == seqID) {
                readthrough += temp.length();
                if (readthrough < start) continue;
                LEN tstart = max(static_cast<LEN>(0), start-(readthrough-static_cast<LEN>(temp.length())));
                if (end > readthrough) {
                    str += temp.substr(tstart);
                } else {
                    str += temp.substr(tstart, static_cast<LEN>(temp.length())-tstart-(readthrough-end));
                    break;
                }
            }
        }
        GetCapitalRNA(str);
        _shift = start;
    }

public:
    /**
     * Sequence length for whole sequence excluding '$'.
     */
    LEN length;
    Sequence() {}
    Sequence(const string& str) : str(str) /* $+AUCU... */ {
        lazy = false;
        part = false;
        _shift = 0;
        length = static_cast<LEN>(str.length())-1;
        SetStrToChar();
    }
    Sequence(const string& str, const vector<char>& sequence, LEN length)
             : str(str), sequence(sequence), length(length) {
                lazy = false;
                part = false;
                _shift = 0;
             }
    Sequence(const string& str, const vector<char>& sequence, LEN length, LEN shift)
             : _shift(shift), str(str), sequence(sequence), length(length) {
        lazy = false;
        part = true;
        _shift = 0;
    }
    Sequence(const string& filename, const int seqID, LEN length) : filename(filename), seqID(seqID), length(length) {
        lazy = true;
        part = true;
        _shift = 0;
    }
    ~Sequence() {}
    /**
     * Set sequence.
     */
    void SetStrToChar() {
        sequence = vector<char>();
        transform(str.begin(), str.end(), back_inserter(sequence), Base());
    }
    /**
     * @return ACGU at pos.
     */
    char strget(LEN pos) const {
        return str[pos-_shift];
    }
    /**
     * @return 0-3 at pos.
     */
    int seqget(LEN pos) const {
        return sequence[pos-_shift];
    }
    /**
     * @return Substring of str that starts at position 'start' of length 'slen'.
     */
    string substr(LEN start, LEN slen) const {
        return str.substr(start-_shift, slen);
    }
    /**
     * @return Substring of str that starts at position 'start'.
     */
    string substr(LEN start) const {
        return str.substr(start-_shift);
    }
    /**
     * Delete a base at 'x'th position on original sequence.
     */
    bool Delete(LEN x) {
        x = x; // slide for $;
        str = str.substr(0, x-_shift)+str.substr(x-_shift+1);
        sequence.erase(sequence.begin()+x-_shift);
        length -= static_cast<LEN>(1);
        return true;
    }
    /**
     * Substitute a base into 'c' at 'x'th position on original sequence.
     * @return whether a base is changed.
     */
    bool Sub(LEN x, char c, char i) {
        x = x; // slide for $;
        if (str[x-_shift] == c) return false;
        str[x-_shift] = c;
        sequence[x-_shift] = i;
        return true;
    }
    /**
     * Insert 'c' before 'x'th position on original sequence.
     */
    bool Insert(LEN x, char c, char i) {
        x = x; // slide for $;
        str.insert(str.begin()+x-_shift, c);
        sequence.insert(sequence.begin()+x-_shift, i);
        length += static_cast<LEN>(1);
        return true;
    }
    /**
     * Cuts out a sequence from 'start' to 'end'.
     */
    void CutSequence(LEN start, LEN end) {
        if (!lazy) {
            if (!part) {
                part = true;
                CutVector(start, end+1);
                CutString(start, end+1);
                _shift = start;
                cout << "#-Cut sequence ----" << start << "-----" << end << "----" << length << endl;
            }
        } else {
            ReadPartialSeq(start, end+1);
            SetStrToChar();
            cout << "#-Read sequence ----" << start << "-----" << end << "----" << length << endl;
        }
        cout << "# " << str.substr(0, 50);
        if (str.length() > 50) cout << "...";
        cout << endl;
    }
    /**
     * @return bp index of base pairing between 'i' and 'j'.
     */
    int slidebp(LEN i, LEN j) const {
        return bp(i-_shift, j-_shift, sequence);
    }
    /**
     * @return Reversed type of base pair between 'i' and 'j'.
     */
    int sliderbp(LEN i, LEN j) const {
        return rtype[bp(i-_shift, j-_shift, sequence)];
    }
    void printSeq() const {
        cout << "#" << _shift << " " << sequence.size() << endl;
        for (LEN i = _shift; i < sequence.size(); i++)
            cout << seqget(i) << " ";
        cout << endl;
    }
};

}
}

#endif