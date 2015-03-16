#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <cmath>
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058
#endif

#ifndef _PLOT_STRUCT_HH
#define _PLOT_STRUCT_HH
#include "matrix.hh"
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#define PROF_NUM (6)

namespace Rfold {

using std::ofstream;
using std::string;
using std::vector;
using std::ostringstream;
using std::next;
using std::fixed;

extern const char* RNAss_head;
extern const char* RNAss_color;
extern const char* RNAss_pcolor;
extern const char *anote_macros;

/**
 * A Plot class.
 * Based on Vienna RNA package 2.0.7 and centroid fold 0.0.10.
 */
class Plot {
private:
	/**
	 * A private variable.
	 * Sequence length.
	 */
	int length;
	/**
	 * A private variable.
	 * Store 0 or an index of base pairing partner.
	 * The first and last index indicates an imaginery base pair.
	 */
	vector<int> pairs;
	/**
	 * A private variable.
	 * Temporary store an angle for each base.
	 */
    Vec angle;
	/**
	 * A private variable.
	 * Store stem probability.
	 */
    Vec stem;
   	/**
	 * A private variable.
	 * Store # N*PROF_NUM of RNA profile values.
	 */
    Vec prof;

	void MakePairMat(string structure);
	int GetCoordinates(Vec& X, Vec& Y);
	int GetNextPair(int count, int p, const vector<int>& pairs, vector<int>& remember);
	void GetLoopAngle(int i, int j, const vector<int>& pairs);
	void BendAngleLoop(int begin, DOUBLE polygon, vector<int>& remember);
	void PlotStruct(ofstream& ofs, string& seq, string& structure);
	void DrawData(ofstream& ofs);
	void SetLength(int tlength) {
		length = tlength;
	}
	void SetStem(Vec& vec) {
		stem = vec;
	}
public:
	Plot() {}
	virtual ~Plot() {}
   	/**
	 * Draws naive structure image.
	 * @param seq base string.
	 * @param structure structure stirng in blacket style.
	 * @param filename filename for image output.
	 * @return 1 if success.
	 */
	static int RNAPlot(string seq, string structure, string filename);
   	/**
	 * Draws structure image with stem probability or RNA profiles.
	 * @param seq base string.
	 * @param structure structure stirng in blacket style.
	 * @param filename filename for image output.
	 * @param stem reference of the vector of stem probability or RNA profile.
	 * @return 1 if success.
	 */
	static int RNAColorPlot(string seq, string structure, string filename, Vec& stem);
};
}
#endif
