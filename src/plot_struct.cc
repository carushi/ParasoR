#include "plot_struct.hh"

namespace Rfold {

/**
 * A postscript template public variable.
 * From PS_dot.c in Vienna RNA package 2.0.7.
 */
const char *RNAss_head =
"%%BeginProlog\n"
"/RNAplot 100 dict def\n"
"RNAplot begin\n"
"/fsize  14 def\n"
"/outlinecolor {0.2 setgray} bind def\n"
"/paircolor    {0.2 setgray} bind def\n"
"/seqcolor     {0   setgray} bind def\n"
"/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/max { 2 copy lt { exch } if pop } bind def\n"
"/arccoords { % i j arccoords\n"
"  % puts optimal x1 y1 x2 y2 coordinates used in bezier curves from i to j\n"
"  % onto the stack\n"
"  dup 3 -1 roll dup 4 -1 roll lt dup dup 5 2 roll {exch} if\n"
"  dup 3 -1 roll dup 3 -1 roll exch sub 1 sub dup\n"
"  4 -2 roll 5 -1 roll {exch} if 4 2 roll\n"
"  sequence length dup 2 div exch 3 1 roll lt \n"
"  {exch 5 -1 roll pop 4 -2 roll exch 4 2 roll}\n"
"  { 4 2 roll 5 -1 roll dup 6 1 roll {exch} if\n"
"    4 -2 roll exch pop dup 3 -1 roll dup 4 1 roll\n"
"    exch add 4 -1 roll dup 5 1 roll sub 1 sub\n"
"    5 -1 roll not {4 -2 roll exch 4 2 roll} if\n"
"  }ifelse\n"
"   % compute the scalingfactor and prepare (1-sf) and sf*r\n"
"  2 mul exch cpr 3 1 roll div dup\n"
"  3 -1 roll mul exch 1 exch sub exch\n"
"   % compute the coordinates\n"
"  3 -1 roll 1 sub coor exch get aload pop % get coord for i\n"
"  4 -1 roll dup 5 1 roll mul 3 -1 roll dup 4 1 roll add exch % calculate y1\n"
"  4 -1 roll dup 5 1 roll mul 3 -1 roll dup 4 1 roll add exch % calculate x1\n"
"  5 -1 roll 1 sub coor exch get aload pop % get coord for j\n"
"  % duplicate j coord\n"
"  dup 3 -1 roll dup 4 1 roll exch 8 2 roll\n"
"  6 -1 roll dup 7 1 roll mul 5 -1 roll dup 6 1 roll add exch % calculate y2\n"
"  6 -1 roll mul 5 -1 roll add exch % calculate x2\n"
"  6 -2 roll % reorder\n"
"} bind def\n"
"/drawoutline {\n"
"  gsave outlinecolor newpath\n"
"  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence\n"
"  currentdict /cutpoint known        % check if cutpoint is defined\n"
"  {coor 0 cutpoint getinterval\n"
"   {aload pop lineto} forall         % draw outline of 1st sequence\n"
"   coor cutpoint 1 add get aload pop\n"
"   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence\n"
"   coor cutpoint 1 add coor length cutpoint 1 add sub getinterval\n"
"   {aload pop lineto} forall}        % draw outline of 2nd sequence\n"
"  {coor {aload pop lineto} forall}   % draw outline as a whole\n"
"  ifelse\n"
"  stroke grestore\n"
"} bind def\n"
"/drawpairs {\n"
"  paircolor\n"
"  0.7 setlinewidth\n"
"  [9 3.01] 9 setdash\n"
"  newpath\n"
"  pairs {aload pop\n"
"      currentdict (cpr) known\n"
"      { exch dup\n"
"        coor  exch 1 sub get aload pop moveto\n"
"        exch arccoords curveto\n"
"      }\n"
"      { coor exch 1 sub get aload pop moveto\n"
"        coor exch 1 sub get aload pop lineto\n"
"      }ifelse\n"
"  } forall\n"
"  stroke\n"
"} bind def\n"
"% draw bases\n"
"/drawbases {\n"
"  [] 0 setdash\n"
"  seqcolor\n"
"  0\n"
"  coor {\n"
"    aload pop moveto\n"
"    dup sequence exch 1 getinterval cshow\n"
"    1 add\n"
"  } forall\n"
"  pop\n"
"} bind def\n\n"
"/init {\n"
"  /Helvetica findfont fsize scalefont setfont\n"
"  1 setlinejoin\n"
"  1 setlinecap\n"
"  0.8 setlinewidth\n"
"  72 216 translate\n"
"  % find the coordinate range\n"
"  /xmax -1000 def /xmin 10000 def\n"
"  /ymax -1000 def /ymin 10000 def\n"
"  coor {\n"
"      aload pop\n"
"      dup ymin lt {dup /ymin exch def} if\n"
"      dup ymax gt {/ymax exch def} {pop} ifelse\n"
"      dup xmin lt {dup /xmin exch def} if\n"
"      dup xmax gt {/xmax exch def} {pop} ifelse\n"
"  } forall\n"
"  /size {xmax xmin sub ymax ymin sub max} bind def\n"
"  72 6 mul size div dup scale\n"
"  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
"  translate\n"
"} bind def\n"
"end\n";

/**
 * A postscript template variable for reliability.
 * From ps_plot.cpp in centroid fold 0.0.10.
 */
const char *anote_macros =
"RNAplot begin\n"
"% extra definitions for standard anotations\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/BLACK { 0 0 0 } def\n"
"/RED   { 1 0 0 } def\n"
"/GREEN { 0 1 0 } def\n"
"/BLUE  { 0 0 1 } def\n"
"/WHITE { 1 1 1 } def\n"
"/LabelFont { % font size LabelFont\n"
"  exch findfont exch fsize mul scalefont setfont\n"
"} bind def\n"
"/Label { % i dx dy (text) Label\n"
"  % write text at base i plus offset dx, dy\n"
"  4 3 roll 1 sub coor exch get aload pop moveto\n"
"  3 1 roll fsize mul exch fsize mul exch rmoveto\n"
"  show\n"
"} bind def\n"
"/cmark { % i cmark   draw circle around base i\n"
"  newpath 1 sub coor exch get aload pop\n"
"  fsize 2 div 0 360 arc stroke\n"
"} bind def\n"
"/gmark { % i j c gmark\n"
"  % draw basepair i,j with c counter examples in gray\n"
"  gsave\n"
"  3 min [0 0.33 0.66 0.9] exch get setgray\n"
"  1 sub dup coor exch get aload pop moveto\n"
"  sequence exch 1 getinterval cshow\n"
"  1 sub dup coor exch get aload pop moveto\n"
"  sequence exch 1 getinterval cshow\n"
"  grestore\n"
"} bind def\n"
"/segmark { % f i j lw r g b segmark\n"
"  % mark segment [i,j] with outline width lw and color rgb\n"
"  % use omark and Fomark instead\n"
"  gsave\n"
"  setrgbcolor setlinewidth\n"
"  newpath\n"
"  1 sub exch 1 sub dup\n"
"  coor exch get aload pop moveto\n"
"  currentdict (cpr) known\n"
"  {\n"
"    3 -1 roll dup 4 1 roll dup\n"
"    {\n"
"      3 1 roll dup 3 -1 roll dup\n"
"      4 1 roll exch 5 2 roll exch\n"
"    }\n"
"    {\n"
"      3 1 roll exch\n"
"    } ifelse\n"
"    1 exch { coor exch get aload pop lineto } for\n"
"    {\n"
"      dup 3 1 roll 1 add exch 1 add arccoords pop pop\n"
"      4 2 roll 5 -1 roll coor exch get aload pop curveto\n"
"    } if\n"
"  }\n"
"  {\n"
"    exch 1 exch {\n"
"      coor exch get aload pop lineto\n"
"    } for\n"
"  } ifelse\n"
"  { closepath fill } if  stroke\n"
"  grestore\n"
"} bind def\n"
"/omark { % i j lw r g b omark\n"
"  % stroke segment [i..j] with linewidth lw, color rgb\n"
"  false 7 1 roll segmark\n"
"} bind def\n"
"/Fomark { % i j r g b Fomark\n"
"  % fill segment [i..j] with color rgb\n"
"  % should precede drawbases\n"
"  1 4 1 roll true 7 1 roll segmark\n"
"} bind def\n"
"/BFmark{ % i j k l r g b BFmark\n"
"  % fill block between pairs (i,j) and (k,l) with color rgb\n"
"  % should precede drawbases\n"
"  gsave\n"
"  setrgbcolor\n"
"  newpath\n"
"  currentdict (cpr) known\n"
"  {\n"
"    dup 1 sub coor exch get aload pop moveto % move to l\n"
"    dup 1 sub 4 -1 roll dup 5 1 roll 1 sub 1 exch\n"
"    { coor exch get aload pop lineto } for % lines from l to j\n"
"    3 -1 roll 4 -1 roll dup 5 1 roll arccoords curveto % curve from j to i\n"
"    exch dup 4 -1 roll 1 sub exch 1 sub 1 exch\n"
"    { coor exch get aload pop lineto } for % lines from i to k\n"
"    exch arccoords curveto% curve from k to l\n"
"  }\n"
"  {  exch 4 3 roll exch 1 sub exch 1 sub dup\n"
"     coor exch get aload pop moveto\n"
"     exch 1 exch { coor exch get aload pop lineto } for\n"
"     exch 1 sub exch 1 sub dup\n"
"     coor exch get aload pop lineto\n"
"     exch 1 exch { coor exch get aload pop lineto } for\n"
"  } ifelse\n"
"    closepath fill stroke\n"
"   grestore\n"
"} bind def\n"
"/hsb {\n"
"  dup 0.3 mul 1 exch sub sethsbcolor\n"
"} bind def\n"
"/colorpair { % i j hue sat colorpair\n"
"  % draw basepair i,j in color\n"
"  % 1 index 0.00 ne {\n"
"  gsave\n"
"  newpath\n"
"  hsb\n"
"  fsize setlinewidth\n"
"  currentdict (cpr) known\n"
"  {\n"
"    exch dup\n"
"    coor  exch 1 sub get aload pop moveto\n"
"    exch arccoords curveto\n"
"  }\n"
"  { 1 sub coor exch get aload pop moveto\n"
"    1 sub coor exch get aload pop lineto\n"
"  } ifelse\n"
"   stroke\n"
"   grestore\n"
"   % } if\n"
"} bind def\n"
 "end\n\n";

/**
 * A postscript template variable for reliability.
 * From RNAss_color, pl_plot.cpp in centroid fold 0.0.10.
 */
const char *RNAss_color =
"/range 0.8 def\n"
"/drawreliability {\n"
"  /Smax 1 def\n"
"  0\n"
"  coor {\n"
"    aload pop\n"
"    S 3 index get\n"
"    Smax div range mul\n"
"    invert {range exch sub} if\n"
"    1 1 sethsbcolor\n"
"    newpath\n"
"    fsize 2 div 0 360 arc\n"
"    fill\n"
"    1 add\n"
"  } forall\n"
"} bind def\n"
"/colorbar { %% xloc yloc colorbar -> []\n"
"  /STR 8 string def\n"
"  gsave\n"
"    xmin xmax add size sub 2 div\n"
"    ymin ymax add size sub 2 div translate\n"
"    size dup scale\n"
"    translate\n"
"    0.015 dup scale\n"
"    /tics 64 def\n"
"    gsave\n"
"      10 tics div 1 scale\n"
"      0 1 tics\n"
"      {\n"
"     dup 0 moveto 0.5 add\n"
"     tics div range mul\n"
"     invert {range exch sub} if\n"
"     1 1 sethsbcolor\n"
"     1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill\n"
"      } for\n"
"    grestore\n"
"    0 setgray\n"
"    -0.1 1.01 moveto (0) gsave 0.1 dup scale show grestore\n"
"    10 1.01 moveto Smax STR cvs\n"
"    gsave 0.1 dup scale dup stringwidth pop -2 div 0 rmoveto show grestore\n"
"  grestore\n"
"} bind def\n";

/**
 * A postscript template variable for profile.
 * Based on RNAss_color, ps_plot.cpp in centroid fold 0.0.10.
 */
const char *RNAss_pcolor =
"/range 0.8 def\n"
"/drawreliability {\n"
"  /Smax 1 def\n"
"  0\n"
"  coor {\n"
"    aload pop\n"
"    S 3 index get\n"
"    /stem exch def\n"
"    O 3 index get\n"
"    /outer exch def\n"
"    H 3 index get\n"
"    /hairpin exch def\n"
"    M 3 index get\n"
"    /multi exch def\n"
"    I 3 index get\n"
"    /inter exch def\n"
// "    B 3 index get\n"
// "    /bulge exch get\n"
"    /y exch def\n"
"    /x exch def\n"
"    newpath\n"
"    1 0 0 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div 0 stem 360 mul arc\n"
"    closepath\n"
"    fill\n"
"    0.2 1 0.2 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div stem 360 mul stem outer add 360 mul arc\n"
"    closepath\n"
"    fill\n"
"    0 0.5 0 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div stem outer add 360 mul stem outer add multi add 360 mul arc\n"
"    closepath\n"
"    fill\n"
"    1 0.2 1 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div stem outer add multi add 360 mul stem outer add multi hairpin add add 360 mul arc\n"
"    closepath\n"
"    fill\n"
"    0.4 0.4 1 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div stem outer add multi hairpin add add 360 mul stem outer add multi hairpin add add inter add 360 mul arc\n"
"    closepath\n"
"    fill\n"
"    1 0.7 0.2 setrgbcolor\n"
"    x y moveto\n"
"    x y fsize 2 div stem outer add multi hairpin add add inter add 360 mul 360 arc\n"
"    closepath\n"
"    fill\n"
"    1 add\n"
"  } forall\n"
"} bind def\n"
"/colorbar { %% xloc yloc colorbar -> []\n"
"  /STR 8 string def\n"
"  gsave\n"
"    xmin xmax add size sub 2 div\n"
"    ymin ymax add size sub 2 div translate\n"
"    size dup scale\n"
"    translate\n"
"    0.015 dup scale\n"
"    /tics 64 def\n"
"    gsave\n"
"      10 tics div 1 scale\n"
"      0 1 tics\n"
"      {\n"
"     dup 0 moveto 0.5 add\n"
"     tics div range mul\n"
"     invert {range exch sub} if\n"
"     1 1 sethsbcolor\n"
"     1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill\n"
"      } for\n"
"    grestore\n"
"    0 setgray\n"
"    -0.1 1.01 moveto (0) gsave 0.1 dup scale show grestore\n"
"    10 1.01 moveto Smax STR cvs\n"
"    gsave 0.1 dup scale dup stringwidth pop -2 div 0 rmoveto show grestore\n"
"  grestore\n"
"} bind def\n";

void Plot::MakePairMat(string structure)
{
    vector<int> stack;
    pairs = vector<int>(length+2, 0);
    pairs[0] = length;
    for (int i = 1; i <= length; i++) {
        switch (structure[i-1]) {
            case '(':
                stack.push_back(i); break;
            case ')':
                if (stack.size() > 0) {
                    pairs[i] = stack[stack.size()-1];
                    pairs[pairs[i]] = i;
                    stack.pop_back();
                }
            case '.':
            default:
                pairs[i] = 0;
                break;
        }
    }
}

int Plot::GetCoordinates(Vec& X, Vec& Y)
{
    DOUBLE INIT_ANGLE=0.;     /* initial bending angle */
    DOUBLE INIT_X = 100.;     /* coordinate of first digit */
    DOUBLE INIT_Y = 100.;     /* see above */
    DOUBLE RADIUS =  15.;
    angle = Vec(length+5);
    GetLoopAngle(0, length+1, pairs);

    DOUBLE alpha = INIT_ANGLE;
    X[0]  = INIT_X;
    Y[0]  = INIT_Y;
    for (int i = 1; i <= length; i++) {
        X[i] = X[i-1]+RADIUS*cos(alpha);
        Y[i] = Y[i-1]+RADIUS*sin(alpha);
        alpha += M_PI-angle[i+1];
    }
    return length;
}

int Plot::GetNextPair(int count, int p, const vector<int>& pairs, vector<int>& remember)
{
    int k = p, l = pairs[p];
    int start_k = k, start_l = l;
    count += 2;
    remember.push_back(k);
    remember.push_back(l);
    int ladder = 0;
    for (; ; k++, l--, ladder++) {
        if (pairs[k] != l)
            break;
    }
    int fill = ladder-2;
    if (ladder >= 2) {
        angle[start_k+1+fill] += M_PI/2.;   /*  Loop entries and    */
        angle[start_l-1-fill] += M_PI/2.;   /*  exits get an        */
        angle[start_k]        += M_PI/2.;   /*  additional PI/2.    */
        angle[start_l]        += M_PI/2.;   /*  Why ? (exercise)    */
        if (ladder > 2) {
            for (; fill >= 1; fill--) {
                angle[start_k+fill] = M_PI;    /*  fill in the angles  */
                angle[start_l-fill] = M_PI;    /*  for the backbone    */
            }
        }
    }
    GetLoopAngle(k, l, pairs);
    return count;
}

void Plot::GetLoopAngle(int i, int j, const vector<int>& pairs)
             /* i, j are the positions AFTER the last pair of a stack; i.e
                i-1 and j+1 are paired. */
{
    int    count = 2;   /* counts the VERTICES of a loop polygon; that's
                           NOT necessarily the number of unpaired bases!
                           Upon entry the loop has already 2 vertices, namely
                           the pair i-1/j+1.  */
    int    bubble = 0; /* bubble counts the unpaired digits in loops */

    vector<int> remember(1, 0);
    for (int p = i, q = j+1; p != q; )
    {
        if (!pairs[p] || (p == 0))
            p++, count++, bubble++;
        else {
            count = GetNextPair(count, p, pairs, remember);
            p = pairs[p]+1;
        }
        if (p == q)
            remember.push_back(q);
    }
    BendAngleLoop(max(i-1, 0), M_PI*(count-2)/(DOUBLE)count, remember); /* bending angle in loop polygon */
}

void Plot::BendAngleLoop(int begin, DOUBLE polygon, vector<int>& remember)
{
    int r = remember.size()-1;
    for (int i = 1; i <= r; i += 2) {
        for (int fill = 0; fill <= remember[i]-begin; fill++)
            angle[begin+fill] += polygon;
        if (i <= r)
            begin = remember[i+1];
    }
}

void Plot::PlotStruct(ofstream& ofs, string& seq, string& structure)
{
    DOUBLE xmin, xmax, ymin, ymax, size;
    Vec X = Vec(length+1), Y = Vec(length+1);
    if (GetCoordinates(X, Y) != length)
        return;
    xmin = *min_element(X.begin(), X.end());
    xmax = *max_element(X.begin(), X.end());
    ymin = *min_element(Y.begin(), Y.end());
    ymax = *max_element(Y.begin(), Y.end());
    size = max((xmax-xmin),(ymax-ymin));
    ofs <<  "%%!PS-Adobe-3.0 EPSF-3.0\n"
        <<  "%%%%Creator: ParasoR\n"
        <<  "%%%%Title: RNA Secondary Structure Plot\n"
        <<  "%%%%BoundingBox: 66 210 518 662\n"
        <<  "%%%%DocumentFonts: Helvetica\n"
        <<  "%%%%Pages: 1\n"
        <<  "%%%%EndComments\n\n"
        <<  "%% to switch off outline pairs of sequence comment or\n"
        <<  "%% delete the appropriate line near the end of the file\n\n";
    ofs << RNAss_head;
    ofs << "%%%%EndProlog\n";
    ofs << "RNAplot begin\n"
        <<  "%% data start here\n";

    /* cut_point */
    int pos = seq.find('&');
    if (pos != std::string::npos) {
        seq[pos] = ' ';
        ofs << "/cutpoint " << pos << " def\n";
    }
    /* sequence */
    ofs << "/sequence (\\\n";
    for (int i = 0; i < length; i+=255) {
        ofs << seq.substr(i, min(255, length-i)) << "\\\n";
    }
    ofs << ") def\n";
    /* coordinates */
    ofs << "/coor [\n";
    for (int i = 0; i < length; i++) {
        ostringstream ss;
        ss.precision(8);
        ss << fixed << "[" << X[i] << " " << Y[i] << "]";
        ofs << ss.str() << endl;
        // fprintf(--, "[%3.8f %3.8f]\n", X[i], Y[i]);
    }
    ofs << "] def\n";
    /* base pairs */
    ofs << "/pairs [\n";
    for (int i = 1; i <= length; i++) {
        if (pairs[i] > i)
            ofs << "[" << i << " " << pairs[i] << "]\n";
    }
    ofs << "] def\n\n";
    ofs << "init\n\n";

}

/* plot ps image of RNA structure from Vienna RNA package*/
int Plot::RNAPlot(string seq, string structure, string filename)
{
    Plot plot;
    plot.SetLength(structure.length());
    ofstream ofs(filename);
    plot.MakePairMat(structure);
    plot.PlotStruct(ofs, seq, structure);
    ofs <<  "%% switch off outline pairs or bases by removing these lines\n"
            "drawoutline\n"
            "drawpairs\n"
            "drawbases\n";
    ofs <<  "%% show it\nshowpage\n"
        <<  "end\n"
        <<  "%%%%EOF\n";
    return 1;
}

int Plot::RNAColorPlot(string seq, string structure, string filename, Vec& stem)
{
    Plot plot;
    plot.SetLength(structure.length());
    plot.SetStem(stem);
    ofstream ofs(filename);
    plot.MakePairMat(structure);
    plot.PlotStruct(ofs, seq, structure);
    if (plot.stem.size() > 0)
        plot.DrawData(ofs);
    /* draw the data */
    ofs <<  "%% switch off outline pairs or bases by removing these lines\n"
            "drawoutline\n"
            "drawpairs\n"
            "drawbases\n";
    ofs <<  "%% show it\nshowpage\n"
        <<  "end\n"
        <<  "%%%%EOF\n";
    return 1; /* success */
}

/* colored by base-pairing probabilities from Centroid fold*/
void Plot::DrawData(ofstream& ofs)
{
    if (stem.size() == PROF_NUM*length) {
        ofs << RNAss_pcolor;
        const char names[6] = {'B', 'O', 'H', 'M', 'S', 'I'};
        for (int i = 0; i < PROF_NUM; i++) {
            ofs << "/" << names[i] << " [\n";
            for (int j = 0; j < length; j++) {
                ostringstream ss;
                ss.precision(8);
                ss << fixed << "      " << stem[j*PROF_NUM+i] << "\n";
                ofs << ss.str();
            }
            ofs << "] def\n\n";
        }
        ofs << "drawreliability\n";
    } else if (stem.size() == length) {
        ofs << RNAss_color;
        ofs << "/S [\n";
        stem.insert(stem.begin(), -1);
        for (int i = 1; i <= length; i++) {
            ostringstream ss;
            ss.precision(8);
            ss << fixed << "      " << stem[i] << "\n";
            ofs << ss.str();
        }
        ofs << "] def\n\n";
        ofs << "/invert true def\n";
        ofs << "drawreliability\n";
        ofs << "0.0 0.0 colorbar\n";
    }
}

}
