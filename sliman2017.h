// This is a header file for the slimanalysis research
// Define everything here
#if !defined SLIMANALYSIS_2017_H
#define SLIMANALYSIS_2017_H
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TMath.h>
#include <initializer_list>
#include <THStack.h>
#include <TStyle.h>
#include <TColor.h>

// Define these for making a centered plot easier
#define REG 1

// Correction to wire-chamber positions
const double wc_xca = -1.21427;
const double wc_yca = -10.0635;

const int NUMCHAN = 8;

const vector<int> TIMESLICES { 5, 6, 7, 8, 9};

// Holds the indeces for the sigma tiles to make things easier
vector<unsigned int> sigmas = {0, 1, 2, 7};
vector<unsigned int> fingers = {3, 4, 5, 6};

// Fiducial region: for each tile, 4 points are needed (8 numbers)
// TB: top, bottom; LR: left, right
// x,y BL; x,y BR; x,y TR; x,y TL
float fiducialX[NUMCHAN][4] = {
#if REG
  {-33,  32,  16, -46}, // EJ-260 CHECKED (TR)
  {-26,  31,  31, -43}, // EJ-260 2P CHECKED (TL)
  {-30,  35,  30, -37}, // EJ-200 CHECKED 
  {-32, -20, -30, -42}, // SCSN-81 finger 1
  {-16,  -4, -15, -27}, // SCSN-81 finger 2
  {  0,  13,   2, -11}, // SCSN-81 finger 3
  { 17,  30,  18,   5}, // SCSN-81 finger 4
  {-38,  34,  22, -49}  // SCSN-81 #sigma CHECKED
#else
  {  -28,    29,    18, -39.5}, // EJ-260 CHECKED 
  {  -29,  26.5,  12.5,   -43}, // EJ-260 2P CHECKED
  {  -30,    30,    22,   -38}, // EJ-200 CHECKED 
  {-33.5, -21.5,   -29,   -41}, // SCSN-81 finger 1 CHECKED
  {  -16,    -4,   -15,   -27}, // SCSN-81 finger 2 CHECKED
  {    0,    12,     2,   -11}, // SCSN-81 finger 3 CHECKED
  {   14,  26.5,    18,     5}, // SCSN-81 finger 4 CHECKED
  {-29.5,  29.5,    20,   -39}  // SCSN-81 #sigma CHECKED
#endif
};

float fiducialY[NUMCHAN][4] = {
#if REG
  {-52, -43, 29, 15}, // EJ-260 CHECKED (TR)
  {-54, -39, 42, 27}, // EJ-260 2P CHECKED (TL)
  {-47, -39, 42, 34}, // EJ-200 CHECKED 
  {-50, -48, 25, 23}, // SCSN-81 finger 1
  {-50, -48, 27, 25}, // SCSN-81 finger 2
  {-50, -48, 30, 28}, // SCSN-81 finger 3
  {-46, -44, 32, 30}, // SCSN-81 finger 4
  {-52, -46, 34, 21}  // SCSN-81 #sigma CHECKED
#else
  {-29.5, -17,   32,   19}, // EJ-260 CHECKED
  {  -41, -30,   39, 27.5}, // EJ-260 2P CHECKED
  {  -48, -42,   36,   30}, // EJ-200 CHECKED
  {  -39, -38, 18.5,   17}, // SCSN-81 finger 1 CHECKED
  {  -41, -39,   29,   27}, // SCSN-81 finger 2 CHECKED
  {  -39, -37, 31.5, 29.5}, // SCSN-81 finger 3 CHECKED
  {  -27, -25,   33,   31}, // SCSN-81 finger 4 CHECKED
  {-35.5, -26,   34,   25}  // SCSN-81 #sigma CHECKED
#endif
};

// FLIPPED THIS TOO!!! TL and TR
float anti_fiducialX[NUMCHAN][4] = {
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80},
  { -80, 80, 80, -80}
};

float anti_fiducialY[NUMCHAN][4] = {
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90},
  { 50, 50, 90, 90}
};

// Identify channels we need to use
struct channel {
  int chan;
  int ieta;
  int idepth;
  string name;
};

vector<channel> channels = 
  {{6,23,5,"EJ_260"},
   {9,23,6,"EJ_260_2P"},
   {12,23,7,"EJ_200"},
   // {3,23,4,"Scint_XS"},
   // {1,23,2,"Scint_XF"},
   {4,24,4,"SCSN_81F1"},
   {7,24,5,"SCSN_81F2"},
   {10,24,6,"SCSN_81F3"},
   {13,24,7,"SCSN_81F4"},
   {5,25,4,"SCSN_81S"}
  };

const string entry[NUMCHAN] = {
  "EJ-260",
  "EJ-260 2P",
  "EJ-200",
  // "Scint-X #sigma",
  // "Scint-X finger",
  "SCSN-81 F1",
  "SCSN-81 F2",
  "SCSN-81 F3",
  "SCSN-81 F4",
  "SCSN-81 #sigma"
};

const int color[NUMCHAN] = {
  kBlack,
  kGreen,
  kBlue,
  // kBlue+1,
  // kGreen+1,
  kRed,
  kRed+1,
  kRed+2,
  kRed+3,
  kViolet
};

const int style[NUMCHAN] = {
  kDotted,
  kSolid,
  kSolid,
  // kSolid,
  // kDashed,
  kDotted,
  kSolid,
  kSolid,
  kSolid,
  kDashed
};

const int NCH = 16;
const int NTS = 10;

float rot_fiducialX[NUMCHAN][4];
float rot_fiducialY[NUMCHAN][4];
float thetas[NUMCHAN];

// Distance from wire-chamber A to plastic tiles, in [mm]
const double z_ex = 7300;

// Change this is you move to lxplus
//const char* slim_dir = "~/TB_Analysis_17/DATA/new_SLIM/";
const char* slim_dir = "/data/users/abelloni/CERN_TB_Jul17/SLIM_2/";

// Some more constants, possibly useful for plot beautification
const int CANVAS_SIZE_X = 500;
const int CANVAS_SIZE_Y = 500;

// Declaration of functions defined in sliman2017.C file
void doAlignmentPlots(bool debug = false, const char* dir = slim_dir);
void doMaps(int flag = 63, bool debug = false, const char* dir = slim_dir);
void doEnergy(int flag = 15, bool debug = false, const char* dir = slim_dir);
void doTimeSlice(bool debug = false, const char* dir = slim_dir);

// Some helper functions

double calc_theta(int channel_num) {
  // 200 => BL                                                  
  if ( channel_num == 2)
    return atan(abs(fiducialY[channel_num][1] - fiducialY[channel_num][0]) /
                abs(fiducialX[channel_num][1] - fiducialX[channel_num][0]));
  // 260, Fingers1, 2, 3, 81S => TR
  else if (channel_num == 0 || channel_num == 3 || channel_num == 4
           || channel_num == 5 || channel_num == 7)
    return atan(abs(fiducialY[channel_num][3] - fiducialY[channel_num][2]) /
                abs(fiducialX[channel_num][3] - fiducialX[channel_num][2]));
  // 260-2p and Finger 4 => TL
  else if (channel_num == 1 || channel_num == 6)
    return atan(abs(fiducialY[channel_num][3] - fiducialY[channel_num][2]) /
                abs(fiducialX[channel_num][3] - fiducialX[channel_num][2]));
  return 0;
}

bool isFiducial(int i, float x_hit, float y_hit) {
  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((fiducialY[i][k]< y_hit && fiducialY[i][j]>=y_hit) ||
         (fiducialY[i][j]< y_hit && fiducialY[i][k]>=y_hit)) &&
        (fiducialX[i][k]<=x_hit || fiducialX[i][j]<=x_hit)) {
      oddNodes^=(fiducialX[i][k]+(y_hit-fiducialY[i][k])/
                 (fiducialY[i][j]-fiducialY[i][k])*
                 (fiducialX[i][j]-fiducialX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}

bool isRotFiducial(int i, float x_hit, float y_hit) {
  if (rot_fiducialX[0][0] == 0 || rot_fiducialY[0][0] == 0)
    cout << "DANGER WILL ROBINSON" << endl;

  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((rot_fiducialY[i][k]< y_hit && rot_fiducialY[i][j]>=y_hit) ||
         (rot_fiducialY[i][j]< y_hit && rot_fiducialY[i][k]>=y_hit)) &&
        (rot_fiducialX[i][k]<=x_hit || rot_fiducialX[i][j]<=x_hit)) {
      oddNodes^=(rot_fiducialX[i][k]+(y_hit-rot_fiducialY[i][k])/
                 (rot_fiducialY[i][j]-rot_fiducialY[i][k])*
                 (rot_fiducialX[i][j]-rot_fiducialX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}

bool isOtherFiducial(int i, float x_hit, float y_hit, float arrX[4][4], float arrY[4][4]) {
  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((arrY[i][k]< y_hit && arrY[i][j]>=y_hit) ||
         (arrY[i][j]< y_hit && arrY[i][k]>=y_hit)) &&
        (arrX[i][k]<=x_hit || arrX[i][j]<=x_hit)) {
      oddNodes^=(arrX[i][k]+(y_hit-arrY[i][k])/
                 (arrY[i][j]-arrY[i][k])*
                 (arrX[i][j]-arrX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}

double rotate_Point(double point_X, double point_Y, int channel_num, char xy) {
  if (toupper(xy) == 'X') {
    return ((point_X*cos(thetas[channel_num])) + (point_Y*sin(thetas[channel_num])));
  }
  else {
    return ((-point_X*sin(thetas[channel_num])) + (point_Y*cos(thetas[channel_num])));
  }
}

bool above_or_below(double px, double py, int chan) {
  double slope = 0, del = 0;
  // Definately Inside
  if (py <= min(rot_fiducialY[chan][2], rot_fiducialY[chan][3]) && 
      py >= max(rot_fiducialY[chan][0], rot_fiducialY[chan][1]))
    return true;
  // Definately Outside
  else if (py > max(rot_fiducialY[chan][2], rot_fiducialY[chan][3]) ||
	   py < min(rot_fiducialY[chan][0], rot_fiducialY[chan][1]))
    return false;
  // Maybe Top Portion
  else if (py >= min( rot_fiducialY[chan][2], rot_fiducialY[chan][3]) &&
	   py <= max( rot_fiducialY[chan][2], rot_fiducialY[chan][3])) {
    slope = (rot_fiducialY[chan][3] - rot_fiducialY[chan][2]) / 
      (rot_fiducialX[chan][3] - rot_fiducialX[chan][2]);
    del = (py - rot_fiducialY[chan][2]) / (px - rot_fiducialX[chan][2]);
    return del <= slope;
  }
  // Maybe Bottom Portion
  else if (py >= min(rot_fiducialY[chan][0], rot_fiducialY[chan][1]) &&
	   py <= max(rot_fiducialY[chan][0], rot_fiducialY[chan][1])){
    slope = (rot_fiducialY[chan][1] - rot_fiducialY[chan][0]) / 
      (rot_fiducialX[chan][1] - rot_fiducialX[chan][0]);
    del = (py - rot_fiducialY[chan][0]) / (px - rot_fiducialX[chan][0]);
    return del >= slope;
  }
  else {
    cout << "ERROR T or B" << endl;
    return false;
  }
}

bool left_or_right( double px, double py, int chan) {
  double slope = 0, del = 0;
  // Definately Inside
  if (px <= min(rot_fiducialX[chan][2], rot_fiducialX[chan][1]) && 
      px >= max(rot_fiducialX[chan][3], rot_fiducialX[chan][0]))
    return true;
  // Definately Outside
  else if (px > max(rot_fiducialX[chan][2], rot_fiducialX[chan][1]) ||
	   px < min(rot_fiducialX[chan][3], rot_fiducialX[chan][0]))
    return false;
  // Maybe Right Portion
  else if (px >= min(rot_fiducialX[chan][2], rot_fiducialX[chan][1]) &&
	   px <= max(rot_fiducialX[chan][2], rot_fiducialX[chan][1])){
    slope = (rot_fiducialY[chan][3] - rot_fiducialY[chan][0]) / 
      (rot_fiducialX[chan][3] - rot_fiducialX[chan][0]);
    del = (py - rot_fiducialY[chan][0]) / (px - rot_fiducialX[chan][0]);
    return del <= slope;
  }
  // Maybe Left Portion
  else if (px >= min(rot_fiducialX[chan][0], rot_fiducialX[chan][3]) &&
           px <= max(rot_fiducialX[chan][0], rot_fiducialX[chan][3])){
    slope = (rot_fiducialY[chan][3] - rot_fiducialY[chan][0]) / 
      (rot_fiducialX[chan][3] - rot_fiducialX[chan][0]);
    del = (py - rot_fiducialY[chan][0]) / (px - rot_fiducialX[chan][0]);
    return del <= slope;
  }
  else {
    cout << "ERROR L or R" << endl;
    return false;
  }
}

// Fills the array thetas with the appropriate angle to use for each tile
void fill_Rot_Array() {
  for (int chan = 0; chan < NUMCHAN; chan++) {
    thetas[chan] = calc_theta(chan); // fill the theta array
    for (int pt = 0; pt < 4; pt++) {
      // Fill the Rotated Arrays
      rot_fiducialX[chan][pt] = rotate_Point(fiducialX[chan][pt], fiducialY[chan][pt]
					     , chan, 'x');
      rot_fiducialY[chan][pt] = rotate_Point(fiducialX[chan][pt], fiducialY[chan][pt]
					     , chan, 'y');
    } // points
  } // channels 
}

// These are the pre-set bins
const double edges[248] = {
  1.58,   4.73,   7.88,   11.0,   14.2,   17.3,   20.5,   23.6,
  26.8,   29.9,   33.1,   36.2,   39.4,   42.5,   45.7,   48.8,
  53.6,   60.1,   66.6,   73.0,   79.5,   86.0,   92.5,   98.9,
  105,    112,    118,    125,    131,    138,    144,    151,
  157,    164,    170,    177,    186,    199,    212,    225,
  238,    251,    264,    277,    289,    302,    315,    328,
  341,    354,    367,    380,    393,    406,    418,    431,
  444,    464,    490,    516,    542,    568,    594,    620,
  645,    670,    695,    720,    745,
  771,    796,    821,    846,    871,    897,    922,    947,
  960,    1010,   1060,   1120,   1170,   1220,   1270,   1320,
  1370,   1430,   1480,   1530,   1580,   1630,   1690,   1740,
  1790,   1840,   1890,   1940,   2020,   2120,   2230,   2330,
  2430,   2540,   2640,   2740,   2850,   2950,   3050,   3150,
  3260,   3360,   3460,   3570,   3670,   3770,   3880,   3980,
  4080,   4240,   4450,   4650,   4860,   5070,   5280,   5490,
  5680,   5880,   6080,   6280,   6480,
  6680,   6890,   7090,   7290,   7490,   7690,   7890,   8090,
  8400,   8810,   9220,   9630,   10000,  10400,  10900,  11300,
  11700,  12100,  12500,  12900,  13300,  13700,  14100,  14500,
  15000,  15400,  15800,  16200,  16800,  17600,  18400,  19300,
  20100,  20900,  21700,  22500,  23400,  24200,  25000,  25800,
  26600,  27500,  28300,  29100,  29900,  30700,  31600,  32400,
  33200,  34400,  36100,  37700,  39400,  41000,  42700,  44300,
  45900,  47600,  49200,  50800,  52500,
  54100,  55700,  57400,  59000,  60600,  62200,  63900,  65500,
  68000,  71300,  74700,  78000,  81400,  84700,  88000,  91400,
  94700,  98100,  101000, 105000, 108000, 111000, 115000, 118000,
  121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,
  168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,
  230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,
  291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000, 398000};


#endif
