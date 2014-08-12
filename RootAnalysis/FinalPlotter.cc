#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>


#include "TFile.h"
#include "TTree.h"

#include "LittlePlotter.h"
#include "TMVA/Reader.h"

void bookPlots(LittlePlotter& plotter);
TMVA::Reader* bookTopVeto(float& f_mW12, float& f_mt12, float& f_mW34, float& f_mt34);
TMVA::Reader* bookMelaMVA(float& f_mX, float& f_yX, float& f_ptX, float& f_m12, float& f_m34, float& abs_cosThetaStar, float& f_cosTheta1, float& f_cosTheta2, float& f_phi, float& modPhi1);
std::string formatNumberForTable(float num);
//void makeTMVAInput(float& f_mX, float& f_yX, float& f_ptX, float& f_m12, float& f_m34, float& abs_cosThetaStar, float& f_cosTheta1, float& f_cosTheta2, float& f_phi, float& modPhi1);
void printCutFlow(LittlePlotter& plotter, std::vector<TString>& categories, std::string caption);
void setupFileList(std::vector<TFile*>& files);
std::map<TString, TTree*> setupOutputTrees(const std::vector<TString>& categories, float& weight,
												float& f_mX, float& f_yX, float& f_ptX, 
												float& f_m12, float& f_m34, float& abs_cosThetaStar, 
												float& f_cosTheta1, float& f_cosTheta2, float& f_phi, float& modPhi1);
void setupTrees(const std::vector<TFile*>& files, std::vector<TString>& categories, std::map<TString, TTree*>& trees, std::map<TString, TTree*>& metaTrees);

//                         l                    c    b                                                tau
double oldTagWeights[] = {0.01, -99, -99, -99, 0.2, 0.7, -99, -99, -99, -99, -99, -99, -99, -99, -99, 0.2};
double newTagWeights[] = {0.01, -99, -99, -99, 0.1, 0.7, -99, -99, -99, -99, -99, -99, -99, -99, -99, 0.2};

int main( int argc, char** argv )
{
	bool makeTmvaInput = false;
	for(int i = 1; i < argc; ++i)
	{
		if(std::string(argv[i]) == "--makeTmvaInput") makeTmvaInput = true;
	}
	std::map<TString, TTree*> outTrees;
	std::vector<TFile*> files;
	setupFileList(files);

	std::vector<TString> categories;	
	std::map<TString, TTree*> trees; //Contains event information
	std::map<TString, TTree*> metaTrees; //Contains information about the jobs used to produce "trees".
	setupTrees(files, categories, trees, metaTrees);

	//Setup variables to pass to Top Veto TMVA
	double ptX, etaX, mX, yX;
	double pt12, m12, mW12, mt12, dRW12;
	double pt34, m34, mW34, mt34, dRW34;
	double cosThetaStar, cosTheta1, cosTheta2, phi, phi1;

	//Annoyingly, TMVA::Reader can only use floats, not doubles! So, duplicate everything.
	float abs_cosThetaStar, abs_cosTheta1, abs_cosTheta2, modPhi1, f_phi;	
	float f_ptX, f_yX, f_mX, f_cosTheta1, f_cosTheta2;
	float f_mW12, f_mt12, f_dRW12, f_m12;
	float f_mW34, f_mt34, f_dRW34, f_m34;
	float weight; 
	TFile* fOut;
	if(makeTmvaInput) 
	{
		fOut = TFile::Open("TMVAInput.root", "RECREATE");
		outTrees = setupOutputTrees(categories, weight, f_mX, f_yX, f_ptX, f_m12, f_m34, abs_cosThetaStar, f_cosTheta1, f_cosTheta2, f_phi, modPhi1);
	}
	//Initialise the TMVA readers for the top veto and the final kinematic selection
	TMVA::Reader* topVeto = bookTopVeto(f_mW12, f_mt12, f_mW34, f_mt34);
	TMVA::Reader* finalMVA = makeTmvaInput ? 0 : bookMelaMVA(f_mX, f_yX, f_ptX, f_m12, f_m34, abs_cosThetaStar, f_cosTheta1, f_cosTheta2, f_phi, modPhi1);

	LittlePlotter plotter(categories);
	bookPlots(plotter);

	for(std::map<TString, TTree*>::iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
	{
		plotter.setCurrentCat(treeIt->first);

		//Establish sample metadata properties before analysing the event tree.
		double sumInputWeights = 0; //This is the total number of events processed, used to calculate xsec weights.
		double sumInputWeightsPerJob;
		TTree* metaTree = metaTrees[treeIt->first];
		long nEntries = metaTree->GetEntries();
		//Don't use nEvent_input, since the Sherpa samples have event weights that != unity
		//metaTree->SetBranchAddress("nEvent_input", nEventsPerJob);
		metaTree->SetBranchAddress("sumW_input", &sumInputWeightsPerJob);
		for(int i = 0; i < nEntries; ++i)
		{
			metaTree->GetEntry(i);
			sumInputWeights += sumInputWeightsPerJob;
		}
		//Calculate weights to account for sample statistics and process cross-sections
		//Cross-sections are in fb.
		float xsec = -99.; 
		if(treeIt->first == "bbbb") xsec = 146030;// * 1.5; 
		else if(treeIt->first == "bbcc") xsec = 317780;// * 1.5; 
		else if(treeIt->first == "ttbar") xsec = 212070;
		else if(treeIt->first == "Hjj") xsec = 3493.9;
		else if(treeIt->first == "Hbb") xsec = 489.24;
		else if(treeIt->first == "ZH") xsec = 35.497;
		else if(treeIt->first == "ttH") xsec = 135.31;
		else if(treeIt->first == "HH") xsec = 11.586;
		//Check we have found a cross-section for this process.	
		if(xsec < -1.)
		{
			std::cout<<"ERROR! Did not find cross-section for process "<< treeIt->first <<", skipping sample."<< std::endl;
			continue;
		}
		double xsecWeight = (xsec*3000.)/sumInputWeights; //N expected events for 3000 fb^{-1}/N events in sample.
		std::cout<< treeIt->first <<": (xsec*3000)/sumInputWeights = ("<< xsec <<"*3000)/"<< sumInputWeights <<" = "<< xsecWeight << std::endl;
		//Now analyse the event tree.
		std::cout<<"Analysing "<< treeIt->first << std::endl;
		TTree* tree = treeIt->second;
		long nEvents = tree->GetEntries();
		double treeWeight, genWeight, btagWeight;
		int flav1, flav2, flav3, flav4;
		tree->SetBranchAddress("weight", &treeWeight);
		tree->SetBranchAddress("genWeight", &genWeight);
		tree->SetBranchAddress("btagWeight", &btagWeight);
		tree->SetBranchAddress("flav1", &flav1);
		tree->SetBranchAddress("flav2", &flav2);
		tree->SetBranchAddress("flav3", &flav3);
		tree->SetBranchAddress("flav4", &flav4);
		tree->SetBranchAddress("yX", &yX);
		tree->SetBranchAddress("ptX", &ptX);
		tree->SetBranchAddress("mX", &mX);
		tree->SetBranchAddress("m12", &m12);
		tree->SetBranchAddress("mW12", &mW12);
		tree->SetBranchAddress("mt12", &mt12);
		tree->SetBranchAddress("dRW12", &dRW12);
		tree->SetBranchAddress("m34", &m34);
		tree->SetBranchAddress("mW34", &mW34);
		tree->SetBranchAddress("mt34", &mt34);
		tree->SetBranchAddress("dRW34", &dRW34);
		tree->SetBranchAddress("pt12", &pt12);
		tree->SetBranchAddress("pt34", &pt34);
		tree->SetBranchAddress("cosThetaStar", &cosThetaStar);
		tree->SetBranchAddress("cosTheta1", &cosTheta1);
		tree->SetBranchAddress("cosTheta2", &cosTheta2);
		tree->SetBranchAddress("Phi", &phi);
		tree->SetBranchAddress("Phi1", &phi1);

		for(int i = 0; i < nEvents; ++i)
		{
            if(i%(nEvents/10) == 0) std::cout<<"Analysing event "<< i <<"/"<< nEvents << std::endl;
			tree->GetEntry(i);

			// Modify the b-tagging weight
			assert(fabs((treeWeight - genWeight*btagWeight)/treeWeight) < 1e-3); // Make sure that we understand the current weight, otherwise abort
			double oldWeight = oldTagWeights[flav1] * oldTagWeights[flav2] * oldTagWeights[flav3] * oldTagWeights[flav4];
			assert(fabs(oldWeight - btagWeight) < 1e-3);
			double newWeight = newTagWeights[flav1] * newTagWeights[flav2] * newTagWeights[flav3] * newTagWeights[flav4];
			assert(newWeight > 0 && newWeight < 1);
			treeWeight *= newWeight/oldWeight;
			//

			weight = xsecWeight * treeWeight;
			abs_cosThetaStar = fabs(cosThetaStar);
			plotter.fill("dijets_absCosThetaStar", abs_cosThetaStar, weight);
			plotter.fill("input_m12", m12, weight);
			plotter.fill("input_m34", m34, weight);
			if(fabs(m12 - 115) > 25. || fabs(m34 - 110.) > 25.) continue;
			plotter.fill("dijets_mH_absCosThetaStar", abs_cosThetaStar, weight);
			f_mW12 = mW12; f_mt12 = mt12; f_dRW12 = dRW12; f_mW34 = mW34; f_mt34 = mt34; f_dRW34 = dRW34;
			//Apply top veto.
			float topMVA = topVeto->EvaluateMVA("BDT");
			plotter.fill("TopVetoBDT", topMVA, weight);
			if(topMVA < -0.2) continue;
			plotter.fill("topVeto_absCosThetaStar", abs_cosThetaStar, weight);
			f_ptX = ptX; f_mX = mX; f_phi = phi; f_yX = yX;
			f_m12 = m12; f_m34 = m34;
			//abs_cosTheta1 = fabs(cosTheta1);
			//abs_cosTheta2 = fabs(cosTheta2);
			f_cosTheta1 = cosTheta1;
			f_cosTheta2 = cosTheta2;
			modPhi1 = fabs(phi1) - 0.5*M_PI;
			if(makeTmvaInput) outTrees[treeIt->first]->Fill();
			plotter.fill("mX", mX, weight);
			plotter.fill("m12", m12, weight);
			plotter.fill("m34", m34, weight);
			plotter.fill("absCosThetaStar", abs_cosThetaStar, weight);
			plotter.fill("cosTheta1", cosTheta1, weight);
			plotter.fill("cosTheta2", cosTheta2, weight);
			plotter.fill("Phi", phi, weight);
			plotter.fill("Phi1", phi1, weight);
			plotter.fill("yX", yX, weight);
			plotter.fill("ptX", ptX, weight);
			float bdt = finalMVA ? finalMVA->EvaluateMVA("BDT") : 0;
			plotter.fill("BDT", bdt, weight);
			if(bdt > 0.16) plotter.fill("postBDT_m34", m34, weight);
			plotter.fill("MelaMVA_absCosThetaStar", abs_cosThetaStar, weight);
		}
	}

	std::vector<TString> signal; signal.push_back("HH");
	std::vector<TString> background; background.push_back("bbbb"); background.push_back("bbcc"); background.push_back("ttbar"); 
	background.push_back("Hjj"); background.push_back("ZH"); background.push_back("ttH"); 

	plotter.plotAlone("BDT", categories);
	plotter.plotAlone("TopVetoBDT", categories);
	plotter.plotAlone("mX", categories);
	plotter.plotAlone("m12", categories);
	plotter.plotAlone("m34", categories);
	plotter.plotAlone("input_m12", categories);
	plotter.plotAlone("dijets_absCosThetaStar", categories);
	plotter.plotAlone("input_m34", categories);
	plotter.plotAlone("absCosThetaStar", categories);
	plotter.plotAlone("cosTheta1", categories);
	plotter.plotAlone("cosTheta2", categories);
	plotter.plotAlone("Phi", categories);
	plotter.plotAlone("Phi1", categories);
	plotter.plotAlone("ptX", categories);
	plotter.plotAlone("yX", categories);
	plotter.plotAlone("postBDT_m34", categories);
	plotter.plotStack("BDT", signal, background);
	plotter.plotSoBVsEff("BDT", signal, background);

	std::vector<TString> qcdBack; qcdBack.push_back("HH"); qcdBack.push_back("bbbb"); qcdBack.push_back("bbcc"); qcdBack.push_back("ttbar"); 
	std::vector<TString> ewkBack; ewkBack.push_back("HH"); ewkBack.push_back("Hjj"); ewkBack.push_back("ZH"); ewkBack.push_back("ttH"); 
	std::cout<<"--------------------Cut flow for QCD backgrounds:-----------------------------------------------------------"<<std::endl;
	printCutFlow(plotter, qcdBack, "QCD backgrounds");	
	std::cout<<"------------------------------------------------------------------------------------------------------------"<<std::endl;
	std::cout<<"--------------------Cut flow for EWK backgrounds:-----------------------------------------------------------"<<std::endl;
	printCutFlow(plotter, ewkBack, "EWK backgrounds");	
	std::cout<<"------------------------------------------------------------------------------------------------------------"<<std::endl;
	std::cout<<"--------------------Cut flow for all backgrounds:-----------------------------------------------------------"<<std::endl;
	printCutFlow(plotter, categories, "All backgrounds");	
	std::cout<<"------------------------------------------------------------------------------------------------------------"<<std::endl;
	//if(makeTmvaInput) fOut->Close();
	if(makeTmvaInput)
	{
		for(std::map<TString, TTree*>::iterator treeIt = outTrees.begin(); treeIt != outTrees.end(); ++ treeIt) treeIt->second->Write();
		fOut->Close();
	}
	return 0;
}
void setupFileList(std::vector<TFile*>& files)
{
	files.push_back(TFile::Open("HH.root", "READ"));
	files.push_back(TFile::Open("bbbb.root", "READ"));
	files.push_back(TFile::Open("bbcc.root", "READ"));
	files.push_back(TFile::Open("ttbar.root", "READ"));
	files.push_back(TFile::Open("Hjj.root", "READ"));
	files.push_back(TFile::Open("ZH.root", "READ"));
	files.push_back(TFile::Open("ttH.root", "READ"));
	std::cout<<"setupFileList: Listed "<< files.size() <<" for processing."<< std::endl;
}
void setupTrees(const std::vector<TFile*>& files, std::vector<TString>& categories, std::map<TString, TTree*>& trees, std::map<TString, TTree*>& metaTrees)
{
	std::cout<<"setupTrees: "<< files.size() <<" files will have trees extracted."<< std::endl;
	for(std::vector<TFile*>::const_iterator file = files.begin(); file != files.end(); ++file)
	{
		TString category = (*file)->GetName();
		category.Remove(category.Last('.'), category.Length());
		categories.push_back(category);
		std::cout<<"Processing "<< (*file)->GetName() <<" as "<< category << std::endl;
		trees.insert(std::pair<TString, TTree*>(category, (TTree*) (*file)->Get("TMVAInput")));
		metaTrees.insert(std::pair<TString, TTree*>(category, (TTree*) (*file)->Get("RunInfo")));
	}
	std::cout<<"setupTrees: Found "<< categories.size() <<" samples, with "<< trees.size() <<" event TTrees and "<< metaTrees.size() <<" job information trees."<< std::endl;
}
void bookPlots(LittlePlotter& plotter)
{
	std::cout<<"bookPlots: booking plots now."<< std::endl;
	plotter.printAllCategories();
	//These are used for printCutFlow, since they are bounded between 0 & 1 so the integration of the plot is trustworthy.
	plotter.book(new TH1F("dijets_absCosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, 0., 1.));
	plotter.book(new TH1F("dijets_mH_absCosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, 0., 1.));
	plotter.book(new TH1F("topVeto_absCosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, 0., 1.));
	plotter.book(new TH1F("MelaMVA_absCosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, 0., 1.));
	plotter.book(new TH1F("input_m12", ";m_{12} [GeV];Number of Events", 50, 0., 250.));
	plotter.book(new TH1F("input_m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
	plotter.book(new TH1F("mX", ";m_{X} [GeV];Number of Events", 50, 250., 750.));
	plotter.book(new TH1F("m12", ";m_{12} [GeV];Number of Events", 50, 0., 250.));
	plotter.book(new TH1F("mW12", ";m_{W,12} [GeV];Number of Events", 50, -100., 150.));
	plotter.book(new TH1F("mt12", ";m_{t,12} [GeV];Number of Events", 50, -100., 400.));
	plotter.book(new TH1F("dRW12", ";#DeltaR_{W,12} ;Number of Events", 50, 0., 2.));
	plotter.book(new TH1F("m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
	plotter.book(new TH1F("mW34", ";m_{W,34} [GeV];Number of Events", 50, -100., 150.));
	plotter.book(new TH1F("mt34", ";m_{t,34} [GeV];Number of Events", 50, -100., 400.));
	plotter.book(new TH1F("dRW34", ";#DeltaR_{W,34} ;Number of Events", 50, 0., 2.));
	plotter.book(new TH1F("absCosThetaStar", ";|cos(#theta^{*})|;Number of Events", 50, 0., 1.));
	plotter.book(new TH1F("cosTheta1", ";cos(#theta^{1});Number of Events", 50, -1., 1.));
	plotter.book(new TH1F("cosTheta2", ";cos(#theta^{2});Number of Events", 50, -1., 1.));
	plotter.book(new TH1F("Phi", ";#Phi;Number of Events", 50, -M_PI, M_PI));
	plotter.book(new TH1F("Phi1", ";#Phi_{1};Number of Events", 50, -M_PI, M_PI));
	plotter.book(new TH1F("yX", ";y_{X};Number of Events", 50, -2.5, 2.5));
	plotter.book(new TH1F("ptX", ";X p_{T} [GeV];Number of Events", 50, 0., 250.));

	plotter.book(new TH1F("TopVetoBDT", ";Top Veto BDT Output;Number of Events", 50, -1., 1.));
	plotter.book(new TH1F("BDT", ";BDT Output;Number of Events", 50, -1., 1.));
	plotter.book(new TH1F("postBDT_m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
	std::cout<<"bookPlots: ending."<< std::endl;
}
TMVA::Reader* bookTopVeto(float& f_mW12, float& f_mt12, float& f_mW34, float& f_mt34)
{	
	TMVA::Reader* topVeto = new TMVA::Reader( "!Color:!Silent" );
	topVeto->AddVariable("mW12", &f_mW12);
	topVeto->AddVariable("mt12", &f_mt12);
	topVeto->AddVariable("mW34", &f_mW34);
	topVeto->AddVariable("mt34", &f_mt34);
	topVeto->BookMVA("BDT", "TopVetoWeights/TMVAClassification_BDT.weights.xml");
	return topVeto;
}
TMVA::Reader* bookMelaMVA(float& f_mX, float& f_yX, float& f_ptX, float& f_m12, float& f_m34, float& abs_cosThetaStar, float& f_cosTheta1, float& f_cosTheta2, float& f_phi, float& modPhi1)
{
	TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );
	reader->AddVariable("mX", &f_mX);
	reader->AddVariable("yX", &f_yX);
	reader->AddVariable("ptX", &f_ptX);
	reader->AddVariable("m12", &f_m12);
	reader->AddVariable("m34", &f_m34);
	reader->AddVariable("absCosThetaStar", &abs_cosThetaStar);
	reader->AddVariable("cosTheta1", &f_cosTheta1);
	reader->AddVariable("cosTheta2", &f_cosTheta2);
	reader->AddVariable("Phi", &f_phi);
	reader->AddVariable("modPhi1", &modPhi1);
   	//factory->AddVariable( "cosTheta1 := abs(cosTheta1)",                "|cos(#theta_{1})|", "", 'D' );
    //factory->AddVariable( "Phi",                "#Phi", "", 'D' );
	reader->BookMVA("BDT", "MelaWeights/TMVAClassification_BDT.weights.xml");
	return reader;
}
void printCutFlow(LittlePlotter& plotter, std::vector<TString>& categories, std::string caption)
{
	std::cout<<"\\begin{table}[t]\\begin{center}\\begin{tabular}{l";
	for(int i = 0; i < categories.size(); ++i) std::cout <<"c";
	std::cout<<"}"<<std::endl;
	std::cout<<"Requirement";
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		std::cout <<" & "<< *cat;
	}
	std::cout<<" \\\\\\hline"<< std::endl;
	TString plotName = "dijets_absCosThetaStar";
	std::cout<<"2 dijets"; 
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
	}
	std::cout<<"\\\\"<< std::endl;
	plotName = "dijets_mH_absCosThetaStar";
	std::cout<<"2 dijets $m_H$"; 
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
	}
	std::cout<<"\\\\"<< std::endl;
	//plotName = "BDT";
	plotName = "topVeto_absCosThetaStar";
	std::cout<<"Top Veto"; 
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
	}
	std::cout<<"\\\\"<< std::endl;
	std::cout<<"MVA"; 
	plotName = "MelaMVA_absCosThetaStar";
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
	}
	std::cout<<"\\\\"<< std::endl;
	std::cout<<"\\hline\\end{tabular}\\caption{"<< caption <<"}\\end{center}\\end{table}"<< std::endl;
	std::cout<<std::endl;
	//Now calculate S/B and S/√B (would be more efficient to calculate integrals once.
	std::cout<<"\\begin{table}[t]\\begin{center}\\begin{tabular}{lcccc}"<<std::endl;
	std::cout<<"Requirement & $S$ & $B$ & $S/B$ & $S/\\sqrt{B}$\\\\\\hline";
	std::cout<<" \\\\\\hline"<< std::endl;
	plotName = "dijets_absCosThetaStar";
	float sig = 0.; float back = 0.;
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		if(*cat == "HH") sig = plotter.getPlot(plotName, *cat)->Integral();
		else back += plotter.getPlot(plotName, *cat)->Integral();
	}
	std::cout<<"2 dijets & "<< formatNumberForTable(sig) <<" & "<< formatNumberForTable(back) <<" & "<< formatNumberForTable(sig/back) <<" & "<< formatNumberForTable(sig/sqrt(back)); 
	std::cout<<"\\\\"<< std::endl;
	plotName = "TopVetoBDT";
	back = 0.;
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		if(*cat == "HH") sig = plotter.getPlot(plotName, *cat)->Integral();
		else back += plotter.getPlot(plotName, *cat)->Integral();
	}
	std::cout<<"2 dijets $m_H$ & "<< formatNumberForTable(sig) <<" & "<< formatNumberForTable(back) <<" & "<< formatNumberForTable(sig/back) <<" & "<< formatNumberForTable(sig/sqrt(back)); 
	std::cout<<"\\\\"<< std::endl;
	//plotName = "BDT";
	plotName = "absCosThetaStar";
	back = 0.;
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		if(*cat == "HH") sig = plotter.getPlot(plotName, *cat)->Integral();
		else back += plotter.getPlot(plotName, *cat)->Integral();
	}
	std::cout<<"Top Veto & "<< formatNumberForTable(sig) <<" & "<< formatNumberForTable(back) <<" & "<< formatNumberForTable(sig/back) <<" & "<< formatNumberForTable(sig/sqrt(back)); 
	std::cout<<"\\\\"<< std::endl;
	plotName = "postBDT_m34";
	back = 0.;
	for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
	{
		if(*cat == "HH") sig = plotter.getPlot(plotName, *cat)->Integral();
		else back += plotter.getPlot(plotName, *cat)->Integral();
	}
	std::cout<<"MVA & "<< formatNumberForTable(sig) <<" & "<< formatNumberForTable(back) <<" & "<< formatNumberForTable(sig/back) <<" & "<< formatNumberForTable(sig/sqrt(back)); 
	std::cout<<"\\\\"<< std::endl;
	std::cout<<"\\hline\\end{tabular}\\caption{"<< caption <<"}\\end{center}\\end{table}"<< std::endl;
}
std::map<TString, TTree*> setupOutputTrees(const std::vector<TString>& categories, float& weight,
												float& f_mX, float& f_yX, float& f_ptX, 
												float& f_m12, float& f_m34, float& abs_cosThetaStar, 
												float& f_cosTheta1, float& f_cosTheta2, float& f_phi, float& modPhi1)
{
	std::map<TString, TTree*> outTrees;
	for(std::vector<TString>::const_iterator catIt = categories.begin(); catIt != categories.end(); ++catIt)
	{
		outTrees.insert(std::pair<TString, TTree*>(*catIt, new TTree(*catIt+"Out", "TMVA compatible tree")));
		outTrees[*catIt]->Branch("weight", &weight, "weight/F");
		outTrees[*catIt]->Branch("mX", &f_mX, "mX/F");
		outTrees[*catIt]->Branch("yX", &f_yX, "yX/F");
		outTrees[*catIt]->Branch("ptX", &f_ptX, "ptX/F");
		outTrees[*catIt]->Branch("m12", &f_m12, "m12/F");
		outTrees[*catIt]->Branch("m34", &f_m34, "m34/F");
		outTrees[*catIt]->Branch("absCosThetaStar", &abs_cosThetaStar, "absCosThetaStar/F");
		outTrees[*catIt]->Branch("Phi", &f_phi, "Phi/F");
		outTrees[*catIt]->Branch("modPhi1", &modPhi1, "modPhi1/F");
		outTrees[*catIt]->Branch("cosTheta1", &f_cosTheta1, "cosTheta1/F");
		outTrees[*catIt]->Branch("cosTheta2", &f_cosTheta2, "cosTheta2/F");
	}
	return outTrees;
}
std::string formatNumberForTable(float num)
{
    float significand = num;
    int exponent = 0;
    if(fabs(num) > 1.)
    {
        while(fabs(significand) > 10)
        {
            significand /= 10.;
            ++exponent;
        }
    }else if(fabs(num) > 0.){
        while(fabs(significand) < 1.)
        {
            significand *= 10.;
            --exponent;
        }
    }else{
        std::cout<<"Input = "<< num <<", output = "<< 0 << std::endl;
        return "0";
    }
    std::stringstream outStr;
    if(exponent > 2) outStr<< std::setprecision(3) <<"$"<< significand <<"\\times10^{"<< exponent <<"}$";
    else outStr << std::setprecision(3) << num;
    return outStr.str();
}
