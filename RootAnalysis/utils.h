#ifndef __UTILS_H
#define __UTILS_H

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TVector.h"
#include <sstream>

float deltaPhi(float,float);
float deltaR(float,float,float,float);
int findCutBin(TH1F*,float,int);
void plot2D(TH2F*, bool, TString);
std::vector<TH1F*> makeProjections(TH2*, TString);
void drawHistVec(std::vector<TH1F*>, TString);
void drawHistVec(std::vector<TH1F*>, std::vector<TH1F*>, TString);
void calcFracInRange(TH1F*,std::vector<TH1F*>);
void drawHisto(TH1*,TString);
void drawHistos(TH1*,TH1*,TString,bool normalise=false,bool showMean=false, bool dynamicYRange=false, bool showRMS=false);
void drawHistos(TH1*,TH1*,TH1*,TString,bool normalise=false,bool showMean=false, bool dynamicYRange=false);
void drawHistos(TH1*,TH1*,TH1*,TH1*,TString,bool normalise=false,bool showMean=false);
void drawHistos(TH1*,TH1*,TH1*,TH1*,TH1*,TH1*,TString,bool normalise=false,bool showMean=false, int colOption=0);
void drawHistos(TH1*,TH1*,TString,float,float);
void drawHistos(TH1*,TH1*,TH1*,TString,float,float);
void drawHistos(TH1*,TH1*,TH1*, TH1*,TString,float,float);
void drawRatio(TH1*,TH1*,const char*,TString,bool normalise=false,bool doFit=false,float minY=-99.,float maxY=-99.);
void drawRatio(TH1*,TH1*,TH1*,const char*,TString,bool normalise=false);
void drawRatio(TH1*,TH1*,TH1*,TH1*,const char*,TString,bool normalise=false);
void drawRatioToA(TH1*,TH1*,TH1*,TH1*,const char*,TString);
void drawProfiles(TProfile*, TProfile*, TString, float, float, TString);
void drawProfiles(TProfile*, TProfile*, TProfile*, TString, float, float, TString, bool drawBCRatio=false);
void drawProfiles(TProfile*, TProfile*, TProfile*, TProfile*, TString, float, float, TString);

void drawProfilesWithSys(TProfile*, TProfile*, TProfile*, TProfile*,TProfile*, TString, float, float, TString, int normaliseToBin=0);
void drawHistosWithSys(TH1*,TH1*,TH1*, TH1*,TH1*,TString,float,float,int normaliseToBin=0);
void drawRatioWithSys(TH1*,TH1*,TH1*,TH1*,TH1*,const char*,TString,int normaliseToBin=0);
void drawProfilesWithSys(TProfile*, TProfile*, TProfile*, TProfile*,TString, float, float, TString, int normaliseToBin=0);
void drawHistosWithSys(TH1*,TH1*,TH1*, TH1*,TString,float,float,int normaliseToBin=0);
void drawRatioWithSys(TH1*,TH1*,TH1*,TH1*,const char*,TString,int normaliseToBin=0);
Double_t effSigma(TH1 * hist);

#endif // __UTILS_H
