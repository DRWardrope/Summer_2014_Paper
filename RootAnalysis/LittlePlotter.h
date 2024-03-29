#ifndef LITTLEPLOTTER_H
#define LITTLEPLOTTER_H
#include <cmath>
#include <iostream>
#include <map>
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "utils.h"

class LittlePlotter
{
	public:
		LittlePlotter(TString outputDir=TString("Plots")): m_outputDir(outputDir) {setupColourPalettes(); m_outputDir = "Plots";};
		LittlePlotter(std::vector<TString>& categories, TString outputDir=TString("Plots")): m_categories(categories), m_outputDir(outputDir) {setupColourPalettes();};
		void book(TH1F*);
		std::pair<float, float> calculateSignificance(TString plotName, TString signal, std::vector<TString>& backgrounds, float minVal, float maxVal);
		void fill(TString name, float value, float weight=1.){ getPlot(name, m_currentCat)->Fill(value, weight);}
		TH1F* getPlot(TString& plotName, TString& category);
		std::vector<TH1F*> getPlots(TString& plotName, std::vector<TString>& categories);
		void plotAlone(TString plotName, std::vector<TString>& categories);
		void plotOverlay(TString plotName, std::vector<TString>& categories);
		void plotStack(TString plotName, std::vector<TString>& signal, std::vector<TString>& background, bool sigAtBottom=false);
		void plotSoBVsEff(TString plotName,std::vector<TString>& signal, std::vector<TString>& background);
		void printAllCategories();
		void printAllPlotNames();
		void scalePlots(TString plotName, std::vector<TString>& categories, std::map<TString,float>& numerator, std::map<TString,float>& denominator);
		void setOutputDir(TString outputDir) { m_outputDir = outputDir;}
		void setCurrentCat(const TString& category){ m_currentCat = category; }
		void setupColourPalettes();
		void writeAll();
		TString legendName(const char*);
	private:
		TString m_outputDir;
		std::vector<TString> m_categories;
		TString m_currentCat;
		std::map<TString, std::vector<TH1F*> > m_th1f;
		std::map<TString, std::vector<TH1F*> >::iterator mIt_th1f;
		int m_sigCol[8], m_backCol[8];
};
#endif
