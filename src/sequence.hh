#ifndef _SEQUENCE_HH
#define _SEQUENCE_HH

#include "pair_mat.hh"
#include "energy_const.hh"

#include <vector>
#include <string>

namespace Rfold {
namespace Parameter {

using std::vector;
using std::string;
using std::min;
using std::max;

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
    LEN _shift; // index to be slided;
    string str; // base string;
    vector<char> sequence;

    void CutVector(LEN start, LEN end) {
        sequence = vector<char>(sequence.begin()+start, sequence.begin()+end);
    }
    void CutString(LEN start, LEN end) {
        cout << "-Cut " << start << " " << end << " " << str.length() << endl;
        str = str.substr(start, end-start);
    }

public:
    /**
     * Sequence length for whole sequence.
     */
    LEN length;
    Sequence() {}
    Sequence(const string& str, const vector<char>& sequence, LEN length)
             : str(str), sequence(sequence), length(length) {
                part = false;
                _shift = 0;
             }
    Sequence(const string& str, const vector<char>& sequence, LEN length, LEN shift)
             : _shift(shift), str(str), sequence(sequence), length(length) {
        part = true;
    }
    ~Sequence() {}
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
        x = x+1; // slide for $;
        str = str.substr(0, x-_shift-1)+str.substr(x-_shift+1);
        sequence.erase(sequence.begin()+x-_shift);
        return true;
    }
    /**
     * Substitute a base into 'c' at 'x'th position on original sequence.
     * @return whether a base is changed.
     */
    bool Sub(LEN x, char c, char i) {
        x = x+1; // slide for $;
        if (str[x-_shift] == c) return false;
        str[x-_shift] = c;
        sequence[x-_shift] = i;
        return true;
    }
    /**
     * Insert 'c' before 'x'th position on original sequence.
     */
    bool Insert(LEN x, char c, char i) {
        x = x+1; // slide for $;
        str.insert(str.begin()+x-_shift, c);
        sequence.insert(sequence.begin()+x-_shift, i);
        return true;
    }
    /**
     * Cuts out a sequence from 'start' to 'end'.
     */
    void CutSequence(LEN start, LEN end) {
        if (!part) {
            part = true;
            CutVector(start, end+1);
            CutString(start, end+1);
            _shift = start;
            cout << "-Cut sequence ----" << start << "-----" << end << "----" << endl;
        }
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
};

}
}

#endif