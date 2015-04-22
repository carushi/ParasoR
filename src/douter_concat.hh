#ifndef _DOUTER_CONCAT_HH
#define _DOUTER_CONCAT_HH

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include "energy_const.hh"


namespace Douter_concat {

using std::string;
using std::vector;
using std::copy;
using std::back_inserter;
using std::istreambuf_iterator;
using std::istream_iterator;
using std::ostream_iterator;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::max;

#define FILELINES (100)

static int binSize(string str) {
	if (str == "long") return 16;
	else if (str == "short") return 4;
	else return 8;
}

static int binSize() {
	return sizeof(DOUBLE);
}

static void WriteHead(ofstream& ofs) {
	int h = 1;
	ofs.write((char*)&h, sizeof(int));
}

static bool CheckHeader(ifstream& ifs, int h, string ifile) {
	if (!ifs) {
		cerr << "No file :" << ifile << endl;
    	return true;
    } else if (h == 1)	return true;
    else if (h > 1000000 || h <= 0) {
        cerr << "File format broken :" << ifile << " header :" << h << endl;
        return true;
	}
	return false;
}

static void MoveFile(string dfile, string sfile)
{
	ifstream ifs(sfile.c_str());
    ofstream ofs(dfile.c_str(), ios::trunc);
    istreambuf_iterator<char> in(ifs), eos;
    ostream_iterator<char> out(ofs);
    copy(in, eos, out);
    ifs.close();
    ofs.close();
    remove(sfile.c_str());
}

static bool ConcatBin(string filename)
{
	string tfile = filename+"_new";
	ifstream ifs(filename.c_str(), ios::binary);
	ofstream ofs;
    int h = 0;
    for (int i = 0; ifs.read((char*)&h, sizeof(int)) && !ifs.eof(); i++) {
    	if (CheckHeader(ifs, h, filename)) return (h != 1);
    	DOUBLE value = 0.;
    	if (i == 0) {
	    	ofs.open(tfile.c_str(), ios::binary | ios::trunc);
	    	WriteHead(ofs);
	    }
    	for (int j = 0; j < h; j++) {
    		ifs.read((char*)&value, sizeof(DOUBLE));
    		ofs.write((char*)&value, sizeof(DOUBLE));
    	}
    }
    ifs.close();
    ofs.close();
    MoveFile(filename, tfile);
    return false;
}

static LEN SeekAllData(string ifile, int h)
{
	LEN pos;
	ifstream ifs(ifile.c_str(), ios::binary);
	ifs.seekg(sizeof(int), std::ios_base::beg);
	for (pos = 0; !ifs.eof(); pos++) {
		ifs.seekg(h*sizeof(DOUBLE), std::ios_base::cur);
		ifs.read((char*)&h, sizeof(int));

	}
	return pos;
}

static void WriteLines(ofstream& ofs, LEN start, LEN end, string ifile)
{
	int h = 0;
	DOUBLE value;
	vector<DOUBLE> outers;
	ifstream ifs(ifile.c_str(), ios::binary);
    ifs.seekg(sizeof(int)*start+sizeof(DOUBLE)*(start), std::ios_base::beg);
    for (LEN i = start; i < end; i++) {
    	ifs.read((char*)&h, sizeof(int));
    	vector<DOUBLE> vec;
    	for (int j = 0; j < h; j++) {
	        ifs.read((char*)&(value), sizeof(DOUBLE));
	        vec.push_back(value);
	    }
       	copy(outers.begin(), outers.end(), back_inserter(vec));
    	outers = vec;
	}
	for (vector<DOUBLE>::iterator it = outers.begin(); it != outers.end(); it++)
		ofs.write((char*)&(*it), sizeof(DOUBLE));
}

static bool ConcatBinReverse(string filename)
{
	string tfile = filename+"_new";
	ifstream ifs(filename.c_str(), ios::binary);
    int h = 0;
    ifs.read((char*)&h, sizeof(int));
    if (CheckHeader(ifs, h, filename)) return (h != 1);
    LEN all = SeekAllData(filename, h);
	ofstream ofs(tfile.c_str(), ios::binary | ios::trunc);
	WriteHead(ofs);
    for (LEN i = all; i >= 0; i -= FILELINES) {
    	WriteLines(ofs, max((LEN)0, i-FILELINES), i, filename);
    }
    ifs.close();
    ofs.close();
    MoveFile(filename, tfile);
    return false;
}

static void GetValues(vector<string>& vec, string str, bool reverse)
{
    istringstream iss(str);
    vector<string> words;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(words));
	for (vector<string>::iterator it = words.begin(); it != words.end();) {
		if (*it == "") it = words.erase(it);
		else	++it;
	}
    if (reverse) {
    	copy(vec.begin(), vec.end(), back_inserter(words));
    	vec = words;
    } else copy(words.begin(), words.end(), back_inserter(vec));
}

static void Concat(string filename, bool reverse = false)
{
	string str;
	ifstream ifs(filename.c_str());
	vector<string> outers;
    while (getline(ifs, str)) {
        if (str.length() == 0) continue;
        GetValues(outers, str, reverse);
    }
    ifs.close();
    ofstream ofs(filename.c_str());
    for (vector<string>::iterator it = outers.begin(); it != outers.end(); it++)
    	ofs << *it << "\n";
    ofs << "\n";
    ofs.close();
}

static void ConcatReverse(string filename) {
	Concat(filename, true);
}


}

#endif