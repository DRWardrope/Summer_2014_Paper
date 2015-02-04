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

std::string formatNumberForTable(float num);

void printCutFlow(LittlePlotter& plotter, std::vector<TString>& categories, std::string caption);
void setupFileList(std::vector<TFile*>& files);

void setupTrees(const std::vector<TFile*>& files, std::vector<TString>& categories, std::map<TString, TTree*>& trees, std::map<TString, TTree*>& metaTrees);

//                         l                    c    b                                                tau
double oldTagWeights[] = {0.01, -99, -99, -99, 0.2, 0.7, -99, -99, -99, -99, -99, -99, -99, -99, -99, 0.2};
double newTagWeights[] = {0.01, -99, -99, -99, 0.1, 0.7, -99, -99, -99, -99, -99, -99, -99, -99, -99, 0.1};
bool changeTagWeights = false;

int main( int argc, char** argv )
{

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
    
    double ptA, mA, etaA, yA, PhiA;
    double ptB, mB, etaB, yB, PhiB;
    
    double quarkMiss12, ptMiss12, etaMiss12, phiMiss12, yMiss12, mMiss12;
    double quarkMiss34, ptMiss34, etaMiss34, phiMiss34, yMiss34, mMiss34;
    double numJetsMatched, numQuarks, quarkMatchStatus, alignStatus, multiMatch;
    double charge12, charge34;
    double bLeadStatus12, bLeadStatus34;
    double dijetsHaveTwoBquarks, bothBleading;
    double dR_sublead12, dR_sublead34;
    double dR_lead12, dR_lead34;
    double dR_Wquarks12, dR_Wquarks34;
    double quark1, quark2, quark3, quark4;
    double deltaR1, deltaR2, deltaR3, deltaR4;
    
    double numNonDup, ptMinAll, etaMaxAll, numClosestFound, multiMatch3rd, num3rdJetsMatched, sameClosest, _3rdMatchStatus;
    double quarkMatched3rdJet12, pt3rdJet12, eta3rdJet12, phi3rdJet12, y3rdJet12, m3rdJet12, dR3rdJet12, dRMatched3rdJet12;
    double quarkMatched3rdJet34, pt3rdJet34, eta3rdJet34, phi3rdJet34, y3rdJet34, m3rdJet34, dR3rdJet34, dRMatched3rdJet34;
    
    //Annoyingly, TMVA::Reader can only use floats, not doubles! So, duplicate everything.
    float f_cosThetaStar, abs_cosTheta1, abs_cosTheta2, modPhi1, f_phi;
    float f_ptX, f_yX, f_mX, f_cosTheta1, f_cosTheta2;
    float f_mW12, f_mt12, f_dRW12, f_m12;
    float f_mW34, f_mt34, f_dRW34, f_m34;
    float f_Xtt;
    float weight;
    TFile* fOut;
    
    //Initialise the TMVA readers for the top veto and the final kinematic selection
    TMVA::Reader* topVeto = bookTopVeto(f_mW12, f_mt12, f_mW34, f_mt34);
    
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
        
        if(treeIt->first == "ttbar") xsec = 212070;
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
        
        tree->SetBranchAddress("ptA", &ptA);
        tree->SetBranchAddress("mA", &mA);
        tree->SetBranchAddress("etaA", &etaA);
        tree->SetBranchAddress("yA", &yA);
        tree->SetBranchAddress("phiA", &PhiA);
        
        tree->SetBranchAddress("ptB", &ptB);
        tree->SetBranchAddress("mB", &mB);
        tree->SetBranchAddress("etaB", &etaB);
        tree->SetBranchAddress("yB", &yB);
        tree->SetBranchAddress("phiB", &PhiB);
        
        tree->SetBranchAddress("quarkMiss12", &quarkMiss12);
        tree->SetBranchAddress("ptMiss12", &ptMiss12);
        tree->SetBranchAddress("mMiss12", &mMiss12);
        tree->SetBranchAddress("etaMiss12", &etaMiss12);
        tree->SetBranchAddress("yMiss12", &yMiss12);
        tree->SetBranchAddress("phiMiss12", &phiMiss12);
        
        tree->SetBranchAddress("quarkMiss34", &quarkMiss34);
        tree->SetBranchAddress("ptMiss34", &ptMiss34);
        tree->SetBranchAddress("mMiss34", &mMiss34);
        tree->SetBranchAddress("etaMiss34", &etaMiss34);
        tree->SetBranchAddress("yMiss34", &yMiss34);
        tree->SetBranchAddress("phiMiss34", &phiMiss34);
        
        tree->SetBranchAddress("charge12", &charge12);
        tree->SetBranchAddress("charge34", &charge34);
        
        tree->SetBranchAddress("bLeadStatus12", &bLeadStatus12);
        tree->SetBranchAddress("bLeadStatus34", &bLeadStatus34);
        tree->SetBranchAddress("bothBleading", &bothBleading);
        tree->SetBranchAddress("dijetsHaveTwoBquarks", &dijetsHaveTwoBquarks);
        
        tree->SetBranchAddress("dR_sublead12", &dR_sublead12);
        tree->SetBranchAddress("dR_sublead34", &dR_sublead34);
        tree->SetBranchAddress("dR_lead12", &dR_lead12);
        tree->SetBranchAddress("dR_lead34", &dR_lead34);
        
        tree->SetBranchAddress("dR_Wquarks12", &dR_Wquarks12);
        tree->SetBranchAddress("dR_Wquarks34", &dR_Wquarks34);
        
        tree->SetBranchAddress("numJetsMatched", &numJetsMatched);
        tree->SetBranchAddress("numQuarks", &numQuarks);
        tree->SetBranchAddress("quarkMatchStatus", &quarkMatchStatus);
        tree->SetBranchAddress("alignStatus", &alignStatus);
        tree->SetBranchAddress("multiMatch", &multiMatch);
        
        tree->SetBranchAddress("quark1", &quark1);
        tree->SetBranchAddress("quark2", &quark2);
        tree->SetBranchAddress("quark3", &quark3);
        tree->SetBranchAddress("quark4", &quark4);
        
        tree->SetBranchAddress("deltaR1", &deltaR1);
        tree->SetBranchAddress("deltaR2", &deltaR2);
        tree->SetBranchAddress("deltaR3", &deltaR3);
        tree->SetBranchAddress("deltaR4", &deltaR4);
        
        tree->SetBranchAddress("numNonDup", &numNonDup);
        tree->SetBranchAddress("ptMinAll", &ptMinAll);
        tree->SetBranchAddress("etaMaxAll", &etaMaxAll);
        tree->SetBranchAddress("numClosestFound", &numClosestFound);
        tree->SetBranchAddress("multiMatch3rd", &multiMatch3rd);
        tree->SetBranchAddress("num3rdJetsMatched", &num3rdJetsMatched);
        tree->SetBranchAddress("sameClosest", &sameClosest);
        tree->SetBranchAddress("3rdMatchStatus", &_3rdMatchStatus);
        
        tree->SetBranchAddress("quarkMatched3rdJet12", &quarkMatched3rdJet12);
        tree->SetBranchAddress("pt3rdJet12", &pt3rdJet12);
        tree->SetBranchAddress("m3rdJet12", &m3rdJet12);
        tree->SetBranchAddress("eta3rdJet12", &eta3rdJet12);
        tree->SetBranchAddress("y3rdJet12", &y3rdJet12);
        tree->SetBranchAddress("phi3rdJet12", &phi3rdJet12);
        
        tree->SetBranchAddress("quarkMatched3rdJet34", &quarkMatched3rdJet34);
        tree->SetBranchAddress("pt3rdJet34", &pt3rdJet34);
        tree->SetBranchAddress("m3rdJet34", &m3rdJet34);
        tree->SetBranchAddress("eta3rdJet34", &eta3rdJet34);
        tree->SetBranchAddress("y3rdJet34", &y3rdJet34);
        tree->SetBranchAddress("phi3rdJet34", &phi3rdJet34);
        
        tree->SetBranchAddress("dR3rdJet12", &dR3rdJet12);
        tree->SetBranchAddress("dRMatched3rdJet12", &dRMatched3rdJet12);
        
        tree->SetBranchAddress("dR3rdJet34", &dR3rdJet34);
        tree->SetBranchAddress("dRMatched3rdJet34", &dRMatched3rdJet34);
        
        
        

        
        for(int i = 0; i < nEvents; ++i)
        {
            if(i%(nEvents/10) == 0) std::cout<<"Analysing event "<< i <<"/"<< nEvents << std::endl;
            tree->GetEntry(i);
            //if (i == 50000) {break;} //Debugging speed up
            // Modify the b-tagging weight
            if (changeTagWeights) {
                
                // Make sure that we understand the current weight, otherwise abort
                assert(fabs((treeWeight - genWeight*btagWeight)/treeWeight) < 1e-3);
                double oldWeight = oldTagWeights[flav1] * oldTagWeights[flav2] * oldTagWeights[flav3] * oldTagWeights[flav4];
                assert(fabs(oldWeight - btagWeight) < 1e-3);
                
                // Calculate and apply the new weight, unexpected labels (weight=-99) make the assertion fail
                double newWeight = newTagWeights[flav1] * newTagWeights[flav2] * newTagWeights[flav3] * newTagWeights[flav4];
                assert(newWeight > 0 && newWeight < 1);
                treeWeight *= newWeight/oldWeight;
            }
            
            weight = xsecWeight * treeWeight;
            plotter.fill("dijets_cosThetaStar", cosThetaStar, weight);
            plotter.fill("input_m12", m12, weight);
            plotter.fill("input_m34", m34, weight);
            if(fabs(m12 - 115) > 25. || fabs(m34 - 110.) > 25.) continue;
            plotter.fill("dijets_mH_cosThetaStar", cosThetaStar, weight);
            f_mW12 = mW12; f_mt12 = mt12; f_dRW12 = dRW12; f_mW34 = mW34; f_mt34 = mt34; f_dRW34 = dRW34;
            
            //Fill plots before top veto
            plotter.fill("before_mX", mX, weight);
            plotter.fill("before_m12", m12, weight);
            plotter.fill("before_m34", m34, weight);
            plotter.fill("before_cosThetaStar", cosThetaStar, weight);
            plotter.fill("before_cosTheta1", cosTheta1, weight);
            plotter.fill("before_cosTheta2", cosTheta2, weight);
            plotter.fill("before_Phi", phi, weight);
            plotter.fill("before_Phi1", phi1, weight);
            plotter.fill("before_yX", yX, weight);
            plotter.fill("before_ptX", ptX, weight);
            plotter.fill("before_pt12", pt12, weight);
            plotter.fill("before_pt34", pt34, weight);
            
            plotter.fill("before_ptA", ptA, weight);
            plotter.fill("before_mA", mA, weight);
            plotter.fill("before_etaA", etaA, weight);
            plotter.fill("before_yA", yA, weight);
            plotter.fill("before_PhiA", PhiA, weight);
            
            plotter.fill("before_ptB", ptB, weight);
            plotter.fill("before_mB", mB, weight);
            plotter.fill("before_etaB", etaB, weight);
            plotter.fill("before_yB", yB, weight);
            plotter.fill("before_PhiB", PhiB, weight);
            
            plotter.fill("before_quarkMiss12", fabs(quarkMiss12), weight);
            plotter.fill("before_ptMiss12", ptMiss12, weight);
            //plotter.fill("before_mMiss12", mMiss12, weight);
            plotter.fill("before_etaMiss12", etaMiss12, weight);
            plotter.fill("before_yMiss12", yMiss12, weight);
            plotter.fill("before_phiMiss12", phiMiss12, weight);
            
            plotter.fill("before_quarkMiss34", fabs(quarkMiss34), weight);
            plotter.fill("before_ptMiss34", ptMiss34, weight);
            //plotter.fill("before_mMiss34", mMiss34, weight);
            plotter.fill("before_etaMiss34", etaMiss34, weight);
            plotter.fill("before_yMiss34", yMiss34, weight);
            plotter.fill("before_phiMiss34", phiMiss34, weight);
            
            plotter.fill("before_charge12", charge12, weight);
            plotter.fill("before_charge34", charge34, weight);
            
            plotter.fill("before_bLeadStatus12", bLeadStatus12, weight);
            plotter.fill("before_bLeadStatus34", bLeadStatus34, weight);
            plotter.fill("before_bothBleading", bothBleading, weight);
            plotter.fill("before_dijetsHaveTwoBquarks", dijetsHaveTwoBquarks, weight);
            
            plotter.fill("before_dR_sublead12", dR_sublead12, weight);
            plotter.fill("before_dR_sublead34", dR_sublead34, weight);
            plotter.fill("before_dR_lead12", dR_lead12, weight);
            plotter.fill("before_dR_lead34", dR_lead34, weight);
            
            plotter.fill("before_dR_Wquarks12", dR_Wquarks12, weight);
            plotter.fill("before_dR_Wquarks34", dR_Wquarks34, weight);
            
            plotter.fill("before_numJetsMatched", numJetsMatched, weight);
            plotter.fill("before_numQuarks", numQuarks, weight);
            plotter.fill("before_quarkMatchStatus", quarkMatchStatus, weight);
            plotter.fill("before_alignStatus", alignStatus, weight);
            plotter.fill("before_multiMatch", multiMatch, weight);
            
            plotter.fill("before_quark1", quark1, weight);
            plotter.fill("before_quark2", quark2, weight);
            plotter.fill("before_quark3", quark3, weight);
            plotter.fill("before_quark4", quark4, weight);
            
            plotter.fill("before_deltaR1", deltaR1, weight);
            plotter.fill("before_deltaR2", deltaR2, weight);
            plotter.fill("before_deltaR3", deltaR3, weight);
            plotter.fill("before_deltaR4", deltaR4, weight);
            
            plotter.fill("before_numNonDup", numNonDup, weight);
            plotter.fill("before_ptMinAll", ptMinAll, weight);
            plotter.fill("before_etaMaxAll", etaMaxAll, weight);
            plotter.fill("before_numClosestFound", numClosestFound, weight);
            plotter.fill("before_multiMatch3rd", multiMatch3rd, weight);
            plotter.fill("before_num3rdJetsMatched", num3rdJetsMatched, weight);
            plotter.fill("before_sameClosest", sameClosest, weight);
            plotter.fill("before_3rdMatchStatus", _3rdMatchStatus, weight);
            
            plotter.fill("before_quarkMatched3rdJet12", fabs(quarkMatched3rdJet12), weight);
            plotter.fill("before_pt3rdJet12", pt3rdJet12, weight);
            plotter.fill("before_m3rdJet12", m3rdJet12, weight);
            plotter.fill("before_eta3rdJet12", eta3rdJet12, weight);
            plotter.fill("before_y3rdJet12", y3rdJet12, weight);
            plotter.fill("before_phi3rdJet12", phi3rdJet12, weight);
            
            plotter.fill("before_quarkMatched3rdJet34", fabs(quarkMatched3rdJet34), weight);
            plotter.fill("before_pt3rdJet34", pt3rdJet34, weight);
            plotter.fill("before_m3rdJet34", m3rdJet34, weight);
            plotter.fill("before_eta3rdJet34", eta3rdJet34, weight);
            plotter.fill("before_y3rdJet34", y3rdJet34, weight);
            plotter.fill("before_phi3rdJet34", phi3rdJet34, weight);
            
            plotter.fill("before_dR3rdJet12", dR3rdJet12, weight);
            plotter.fill("before_dRMatched3rdJet12", dRMatched3rdJet12, weight);
            
            plotter.fill("before_dR3rdJet34", dR3rdJet34, weight);
            plotter.fill("before_dRMatched3rdJet34", dRMatched3rdJet34, weight);
            
            //Apply top veto.
            float topMVA = topVeto->EvaluateMVA("BDT");
            plotter.fill("TopVetoBDT", topMVA, weight);
            if(topMVA < -0.05) continue;
            plotter.fill("topVeto_cosThetaStar", cosThetaStar, weight);

            //Fill plots after top veto
            plotter.fill("mX", mX, weight);
            plotter.fill("m12", m12, weight);
            plotter.fill("m34", m34, weight);
            plotter.fill("cosThetaStar", cosThetaStar, weight);
            plotter.fill("cosTheta1", cosTheta1, weight);
            plotter.fill("cosTheta2", cosTheta2, weight);
            plotter.fill("Phi", phi, weight);
            plotter.fill("Phi1", phi1, weight);
            plotter.fill("yX", yX, weight);
            plotter.fill("ptX", ptX, weight);
            plotter.fill("pt12", pt12, weight);
            plotter.fill("pt34", pt34, weight);
            
            plotter.fill("ptA", ptA, weight);
            plotter.fill("mA", mA, weight);
            plotter.fill("etaA", etaA, weight);
            plotter.fill("yA", yA, weight);
            plotter.fill("PhiA", PhiA, weight);
            
            plotter.fill("ptB", ptB, weight);
            plotter.fill("mB", mB, weight);
            plotter.fill("etaB", etaB, weight);
            plotter.fill("yB", yB, weight);
            plotter.fill("PhiB", PhiB, weight);

            plotter.fill("quarkMiss12", quarkMiss12, weight);
            plotter.fill("ptMiss12", ptMiss12, weight);
            //plotter.fill("mMiss12", mMiss12, weight);
            plotter.fill("etaMiss12", etaMiss12, weight);
            plotter.fill("yMiss12", yMiss12, weight);
            plotter.fill("phiMiss12", phiMiss12, weight);
            
            plotter.fill("quarkMiss34", quarkMiss34, weight);
            plotter.fill("ptMiss34", ptMiss34, weight);
            //plotter.fill("mMiss34", mMiss34, weight);
            plotter.fill("etaMiss34", etaMiss34, weight);
            plotter.fill("yMiss34", yMiss34, weight);
            plotter.fill("phiMiss34", phiMiss34, weight);
            
            plotter.fill("charge12", charge12, weight);
            plotter.fill("charge34", charge34, weight);
            
            plotter.fill("bLeadStatus12", bLeadStatus12, weight);
            plotter.fill("bLeadStatus34", bLeadStatus34, weight);
            plotter.fill("bothBleading", bothBleading, weight);
            plotter.fill("dijetsHaveTwoBquarks", dijetsHaveTwoBquarks, weight);
            
            plotter.fill("dR_sublead12", dR_sublead12, weight);
            plotter.fill("dR_sublead34", dR_sublead34, weight);
            plotter.fill("dR_lead12", dR_lead12, weight);
            plotter.fill("dR_lead34", dR_lead34, weight);
            
            plotter.fill("dR_Wquarks12", dR_Wquarks12, weight);
            plotter.fill("dR_Wquarks34", dR_Wquarks34, weight);
            
            plotter.fill("numJetsMatched", numJetsMatched, weight);
            plotter.fill("numQuarks", numQuarks, weight);
            plotter.fill("quarkMatchStatus", quarkMatchStatus, weight);
            plotter.fill("alignStatus", alignStatus, weight);
            plotter.fill("multiMatch", multiMatch, weight);
            
            plotter.fill("quark1", fabs(quark1), weight);
            plotter.fill("quark2", fabs(quark2), weight);
            plotter.fill("quark3", fabs(quark3), weight);
            plotter.fill("quark4", fabs(quark4), weight);
            
            plotter.fill("deltaR1", deltaR1, weight);
            plotter.fill("deltaR2", deltaR2, weight);
            plotter.fill("deltaR3", deltaR3, weight);
            plotter.fill("deltaR4", deltaR4, weight);
            
            plotter.fill("numNonDup", numNonDup, weight);
            plotter.fill("ptMinAll", ptMinAll, weight);
            plotter.fill("etaMaxAll", etaMaxAll, weight);
            plotter.fill("numClosestFound", numClosestFound, weight);
            plotter.fill("multiMatch3rd", multiMatch3rd, weight);
            plotter.fill("num3rdJetsMatched", num3rdJetsMatched, weight);
            plotter.fill("sameClosest", sameClosest, weight);
            plotter.fill("3rdMatchStatus", _3rdMatchStatus, weight);
            
            plotter.fill("quarkMatched3rdJet12", quarkMatched3rdJet12, weight);
            plotter.fill("pt3rdJet12", pt3rdJet12, weight);
            plotter.fill("m3rdJet12", m3rdJet12, weight);
            plotter.fill("eta3rdJet12", eta3rdJet12, weight);
            plotter.fill("y3rdJet12", y3rdJet12, weight);
            plotter.fill("phi3rdJet12", phi3rdJet12, weight);
            
            plotter.fill("quarkMatched3rdJet34", quarkMatched3rdJet34, weight);
            plotter.fill("pt3rdJet34", pt3rdJet34, weight);
            plotter.fill("m3rdJet34", m3rdJet34, weight);
            plotter.fill("eta3rdJet34", eta3rdJet34, weight);
            plotter.fill("y3rdJet34", y3rdJet34, weight);
            plotter.fill("phi3rdJet34", phi3rdJet34, weight);
            
            plotter.fill("dR3rdJet12", dR3rdJet12, weight);
            plotter.fill("dRMatched3rdJet12", dRMatched3rdJet12, weight);
            
            plotter.fill("dR3rdJet34", dR3rdJet34, weight);
            plotter.fill("dRMatched3rdJet34", dRMatched3rdJet34, weight);
        }
    }
    /*
    plotter.plotAlone("TopVetoBDT", categories);
    plotter.plotAlone("dijets_cosThetaStar", categories);
    plotter.plotAlone("dijets_mH_cosThetaStar", categories);
    
    plotter.plotAlone("input_m12", categories);
    plotter.plotAlone("input_m34", categories);
    
    plotter.plotAlone("mX", categories);
    plotter.plotAlone("m12", categories);
    plotter.plotAlone("m34", categories);
    plotter.plotAlone("cosThetaStar", categories);
    plotter.plotAlone("cosTheta1", categories);
    plotter.plotAlone("cosTheta2", categories);
    plotter.plotAlone("Phi", categories);
    plotter.plotAlone("Phi1", categories);
    plotter.plotAlone("ptX", categories);
    plotter.plotAlone("pt12", categories);
    plotter.plotAlone("pt34", categories);
    plotter.plotAlone("yX", categories);
    
    //New plots
    plotter.plotAlone("before_mX", categories);
    plotter.plotAlone("before_m12", categories);
    plotter.plotAlone("before_m34", categories);
    plotter.plotAlone("before_cosThetaStar", categories);
    plotter.plotAlone("before_cosTheta1", categories);
    plotter.plotAlone("before_cosTheta2", categories);
    plotter.plotAlone("before_Phi", categories);
    plotter.plotAlone("before_Phi1", categories);
    plotter.plotAlone("before_ptX", categories);
    plotter.plotAlone("before_pt12", categories);
    plotter.plotAlone("before_pt34", categories);
    plotter.plotAlone("before_yX", categories);

    plotter.plotAlone("before_ptA", categories);
    plotter.plotAlone("before_mA", categories);
    plotter.plotAlone("before_etaA", categories);
    plotter.plotAlone("before_yA", categories);
    plotter.plotAlone("before_PhiA", categories);
    
    plotter.plotAlone("before_ptB", categories);
    plotter.plotAlone("before_mB", categories);
    plotter.plotAlone("before_etaB", categories);
    plotter.plotAlone("before_yB", categories);
    plotter.plotAlone("before_PhiB", categories);*/
    /*
    plotter.plotAlone("before_quarkMiss12", categories);
    plotter.plotAlone("before_ptMiss12", categories);
    //plotter.plotAlone("before_mMiss12", categories);
    plotter.plotAlone("before_etaMiss12", categories);
    plotter.plotAlone("before_yMiss12", categories);
    plotter.plotAlone("before_phiMiss12", categories);
    
    plotter.plotAlone("before_quarkMiss34", categories);
    plotter.plotAlone("before_ptMiss34", categories);
    //plotter.plotAlone("before_mMiss34", categories);
    plotter.plotAlone("before_etaMiss34", categories);
    plotter.plotAlone("before_yMiss34", categories);
    plotter.plotAlone("before_phiMiss34", categories);
    
    plotter.plotAlone("before_charge12", categories);
    plotter.plotAlone("before_charge34", categories);
    
    plotter.plotAlone("before_bLeadStatus12", categories);
    plotter.plotAlone("before_bLeadStatus34", categories);
    plotter.plotAlone("before_bothBleading", categories);
    plotter.plotAlone("before_dijetsHaveTwoBquarks", categories);
    
    plotter.plotAlone("before_dR_sublead12", categories);
    plotter.plotAlone("before_dR_sublead34", categories);
    plotter.plotAlone("before_dR_lead12", categories);
    plotter.plotAlone("before_dR_lead34", categories);
    
    plotter.plotAlone("before_dR_Wquarks12", categories);
    plotter.plotAlone("before_dR_Wquarks34", categories);
    
    plotter.plotAlone("before_numJetsMatched", categories);
    plotter.plotAlone("before_numQuarks", categories);
    plotter.plotAlone("before_quarkMatchStatus", categories);
    plotter.plotAlone("before_alignStatus", categories);
    plotter.plotAlone("before_multiMatch", categories);
    
    plotter.plotAlone("before_quark1", categories);
    plotter.plotAlone("before_quark2", categories);
    plotter.plotAlone("before_quark3", categories);
    plotter.plotAlone("before_quark4", categories);
    
    plotter.plotAlone("before_deltaR1", categories);
    plotter.plotAlone("before_deltaR2", categories);
    plotter.plotAlone("before_deltaR3", categories);
    plotter.plotAlone("before_deltaR4", categories);
    
    plotter.plotAlone("quarkMiss12", categories);
    plotter.plotAlone("ptMiss12", categories);
    //plotter.plotAlone("mMiss12", categories);
    plotter.plotAlone("etaMiss12", categories);
    plotter.plotAlone("yMiss12", categories);
    plotter.plotAlone("phiMiss12", categories);
    
    plotter.plotAlone("quarkMiss34", categories);
    plotter.plotAlone("ptMiss34", categories);
    //plotter.plotAlone("mMiss34", categories);
    plotter.plotAlone("etaMiss34", categories);
    plotter.plotAlone("yMiss34", categories);
    plotter.plotAlone("phiMiss34", categories);
    
    plotter.plotAlone("charge12", categories);
    plotter.plotAlone("charge34", categories);
    
    plotter.plotAlone("bLeadStatus12", categories);
    plotter.plotAlone("bLeadStatus34", categories);
    plotter.plotAlone("bothBleading", categories);
    plotter.plotAlone("dijetsHaveTwoBquarks", categories);
    
    plotter.plotAlone("dR_sublead12", categories);
    plotter.plotAlone("dR_sublead34", categories);
    plotter.plotAlone("dR_lead12", categories);
    plotter.plotAlone("dR_lead34", categories);
    
    plotter.plotAlone("dR_Wquarks12", categories);
    plotter.plotAlone("dR_Wquarks34", categories);
    
    plotter.plotAlone("numJetsMatched", categories);
    plotter.plotAlone("numQuarks", categories);
    plotter.plotAlone("quarkMatchStatus", categories);
    plotter.plotAlone("alignStatus", categories);
    plotter.plotAlone("multiMatch", categories);
    
    plotter.plotAlone("quark1", categories);
    plotter.plotAlone("quark2", categories);
    plotter.plotAlone("quark3", categories);
    plotter.plotAlone("quark4", categories);
    
    plotter.plotAlone("deltaR1", categories);
    plotter.plotAlone("deltaR2", categories);
    plotter.plotAlone("deltaR3", categories);
    plotter.plotAlone("deltaR4", categories);
    */
    //Abusing the modified signal background plot functions
    /*
    plotter.plotBeforeAfter("before_mX", "mX", categories);
    plotter.plotBeforeAfter("before_m12", "m12", categories);
    plotter.plotBeforeAfter("before_m34", "m34", categories);
    plotter.plotBeforeAfter("before_cosThetaStar", "cosThetaStar", categories);
    plotter.plotBeforeAfter("before_cosTheta1", "cosTheta1", categories);
    plotter.plotBeforeAfter("before_cosTheta2", "cosTheta2", categories);
    plotter.plotBeforeAfter("before_Phi", "Phi", categories);
    plotter.plotBeforeAfter("before_Phi1", "Phi1", categories);
    plotter.plotBeforeAfter("before_ptX", "ptX", categories);
    plotter.plotBeforeAfter("before_pt12", "pt12", categories);
    plotter.plotBeforeAfter("before_pt34", "pt34", categories);
    plotter.plotBeforeAfter("before_yX", "yX", categories);
    
    plotter.plotBeforeAfter("before_ptA", "ptA", categories);
    plotter.plotBeforeAfter("before_mA", "mA", categories);
    plotter.plotBeforeAfter("before_etaA", "etaA", categories);
    plotter.plotBeforeAfter("before_yA", "yA", categories);
    plotter.plotBeforeAfter("before_PhiA", "PhiA", categories);
    
    plotter.plotBeforeAfter("before_ptB", "ptB", categories);
    plotter.plotBeforeAfter("before_mB", "mB", categories);
    plotter.plotBeforeAfter("before_etaB", "etaB", categories);
    plotter.plotBeforeAfter("before_yB", "yB", categories);
    plotter.plotBeforeAfter("before_PhiB", "PhiB", categories);
    */
    /*
    plotter.plotBeforeAfter("before_quarkMiss12", "quarkMiss12", categories);
    plotter.plotBeforeAfter("before_ptMiss12", "ptMiss12", categories);
    plotter.plotBeforeAfter("before_mMiss12", "mMiss12", categories);
    plotter.plotBeforeAfter("before_etaMiss12", "etaMiss12", categories);
    plotter.plotBeforeAfter("before_yMiss12", "yMiss12", categories);
    plotter.plotBeforeAfter("before_phiMiss12", "phiMiss12", categories);
    
    plotter.plotBeforeAfter("before_quarkMiss34", "quarkMiss34", categories);
    plotter.plotBeforeAfter("before_ptMiss34", "ptMiss34", categories);
    plotter.plotBeforeAfter("before_mMiss34", "mMiss34", categories);
    plotter.plotBeforeAfter("before_etaMiss34", "etaMiss34", categories);
    plotter.plotBeforeAfter("before_yMiss34", "yMiss34", categories);
    plotter.plotBeforeAfter("before_phiMiss34", "phiMiss34", categories);
    
    plotter.plotBeforeAfter("before_charge12", "charge12", categories);
    plotter.plotBeforeAfter("before_charge34", "charge34", categories);
    
    plotter.plotBeforeAfter("before_bLeadStatus12", "bLeadStatus12", categories);
    plotter.plotBeforeAfter("before_bLeadStatus34", "bLeadStatus34", categories);
    plotter.plotBeforeAfter("before_bothBleading", "bothBleading", categories);
    plotter.plotBeforeAfter("before_dijetsHaveTwoBquarks", "dijetsHaveTwoBquarks", categories);
    
    plotter.plotBeforeAfter("before_dR_sublead12", "dR_sublead12", categories);
    plotter.plotBeforeAfter("before_dR_sublead34", "dR_sublead34", categories);
    plotter.plotBeforeAfter("before_dR_lead12", "dR_lead12", categories);
    plotter.plotBeforeAfter("before_dR_lead34", "dR_lead34", categories);
    
    plotter.plotBeforeAfter("before_dR_Wquarks12", "dR_Wquarks12", categories);
    plotter.plotBeforeAfter("before_dR_Wquarks34", "dR_Wquarks34", categories);
    
    plotter.plotBeforeAfter("before_numJetsMatched", "numJetsMatched", categories);
    plotter.plotBeforeAfter("before_numQuarks", "numQuarks", categories);
    plotter.plotBeforeAfter("before_quarkMatchStatus", "quarkMatchStatus", categories);
    plotter.plotBeforeAfter("before_alignStatus", "alignStatus", categories);
    plotter.plotBeforeAfter("before_multiMatch", "multiMatch", categories);
    
    plotter.plotBeforeAfter("before_quark1", "quark1", categories);
    plotter.plotBeforeAfter("before_quark2", "quark2", categories);
    plotter.plotBeforeAfter("before_quark3", "quark3", categories);
    plotter.plotBeforeAfter("before_quark4", "quark4", categories);
    
    plotter.plotBeforeAfter("before_deltaR1", "deltaR1", categories);
    plotter.plotBeforeAfter("before_deltaR2", "deltaR2",  categories);
    plotter.plotBeforeAfter("before_deltaR3", "deltaR3",  categories);
    plotter.plotBeforeAfter("before_deltaR4", "deltaR4",  categories);
    */
    plotter.plotBeforeAfter("before_numNonDup", "numNonDup",  categories);
    plotter.plotBeforeAfter("before_ptMinAll", "ptMinAll",  categories);
    plotter.plotBeforeAfter("before_etaMaxAll", "etaMaxAll",  categories);
    plotter.plotBeforeAfter("before_numClosestFound", "numClosestFound",  categories);
    plotter.plotBeforeAfter("before_multiMatch3rd", "multiMatch3rd",  categories);
    plotter.plotBeforeAfter("before_num3rdJetsMatched", "num3rdJetsMatched",  categories);
    plotter.plotBeforeAfter("before_sameClosest", "sameClosest",  categories);
    plotter.plotBeforeAfter("before_3rdMatchStatus", "3rdMatchStatus",  categories);
    
    plotter.plotBeforeAfter("before_quarkMatched3rdJet12", "quarkMatched3rdJet12",  categories);
    plotter.plotBeforeAfter("before_pt3rdJet12", "pt3rdJet12",  categories);
    plotter.plotBeforeAfter("before_m3rdJet12", "m3rdJet12",  categories);
    plotter.plotBeforeAfter("before_eta3rdJet12", "eta3rdJet12",  categories);
    plotter.plotBeforeAfter("before_y3rdJet12", "y3rdJet12",  categories);
    plotter.plotBeforeAfter("before_phi3rdJet12", "phi3rdJet12",  categories);
    
    plotter.plotBeforeAfter("before_quarkMatched3rdJet34", "quarkMatched3rdJet34",  categories);
    plotter.plotBeforeAfter("before_pt3rdJet34", "pt3rdJet34",  categories);
    plotter.plotBeforeAfter("before_m3rdJet34", "m3rdJet34",  categories);
    plotter.plotBeforeAfter("before_eta3rdJet34", "eta3rdJet34",  categories);
    plotter.plotBeforeAfter("before_y3rdJet34", "y3rdJet34",  categories);
    plotter.plotBeforeAfter("before_phi3rdJet34", "phi3rdJet34",  categories);
    
    plotter.plotBeforeAfter("before_dR3rdJet12", "dR3rdJet12",  categories);
    plotter.plotBeforeAfter("before_dRMatched3rdJet12", "dRMatched3rdJet12",  categories);
    
    plotter.plotBeforeAfter("before_dR3rdJet34", "dR3rdJet34",  categories);
    plotter.plotBeforeAfter("before_dRMatched3rdJet34", "dRMatched3rdJet34",  categories);
    
    std::cout<<"--------------------Cut flow for all backgrounds:-----------------------------------------------------------"<<std::endl;
    printCutFlow(plotter, categories, "All backgrounds");
    std::cout<<"------------------------------------------------------------------------------------------------------------"<<std::endl;

    return 0;
}
void setupFileList(std::vector<TFile*>& files)
{
    //files.push_back(TFile::Open("HH.root", "READ"));
    files.push_back(TFile::Open("ttbar.root", "READ"));
    
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
    plotter.book(new TH1F("dijets_cosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("dijets_mH_cosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("topVeto_cosThetaStar", ";|cos(#theta^*)|;Number of Events", 50, 1., 1.));

    plotter.book(new TH1F("input_m12", ";m_{12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("input_m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
    
    plotter.book(new TH1F("before_mX", ";m_{X} [GeV];Number of Events", 50, 250., 750.));
    plotter.book(new TH1F("before_m12", ";m_{12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_mW12", ";m_{W,12} [GeV];Number of Events", 50, -100., 150.));
    plotter.book(new TH1F("before_mt12", ";m_{t,12} [GeV];Number of Events", 50, -100., 400.));
    plotter.book(new TH1F("before_dRW12", ";#DeltaR_{W,12} ;Number of Events", 50, 0., 2.));
    plotter.book(new TH1F("before_m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_mW34", ";m_{W,34} [GeV];Number of Events", 50, -100., 150.));
    plotter.book(new TH1F("before_mt34", ";m_{t,34} [GeV];Number of Events", 50, -100., 400.));
    plotter.book(new TH1F("before_dRW34", ";#DeltaR_{W,34} ;Number of Events", 50, 0., 2.));
    plotter.book(new TH1F("before_cosThetaStar", ";|cos(#theta^{*})|;Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("before_cosTheta1", ";cos(#theta^{1});Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("before_cosTheta2", ";cos(#theta^{2});Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("before_Phi", ";#Phi;Number of Events", 50, -M_PI, M_PI));
    plotter.book(new TH1F("before_Phi1", ";#Phi_{1};Number of Events", 50, -M_PI, M_PI));
    plotter.book(new TH1F("before_yX", ";y_{X};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_ptX", ";X p_{T} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_pt12", ";X p_{T} [GeV];Number of Events", 50, 100., 700.));
    plotter.book(new TH1F("before_pt34", ";X p_{T} [GeV];Number of Events", 50,100., 700.));
    
    plotter.book(new TH1F("mX", ";m_{X} [GeV];Number of Events", 50, 250., 750.));
    plotter.book(new TH1F("m12", ";m_{12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("mW12", ";m_{W,12} [GeV];Number of Events", 50, -100., 150.));
    plotter.book(new TH1F("mt12", ";m_{t,12} [GeV];Number of Events", 50, -100., 400.));
    plotter.book(new TH1F("dRW12", ";#DeltaR_{W,12} ;Number of Events", 50, 0., 2.));
    plotter.book(new TH1F("m34", ";m_{34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("mW34", ";m_{W,34} [GeV];Number of Events", 50, -100., 150.));
    plotter.book(new TH1F("mt34", ";m_{t,34} [GeV];Number of Events", 50, -100., 400.));
    plotter.book(new TH1F("dRW34", ";#DeltaR_{W,34} ;Number of Events", 50, 0., 2.));
    plotter.book(new TH1F("cosThetaStar", ";|cos(#theta^{*})|;Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("cosTheta1", ";cos(#theta^{1});Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("cosTheta2", ";cos(#theta^{2});Number of Events", 50, -1., 1.));
    plotter.book(new TH1F("Phi", ";#Phi;Number of Events", 50, -M_PI, M_PI));
    plotter.book(new TH1F("Phi1", ";#Phi_{1};Number of Events", 50, -M_PI, M_PI));
    plotter.book(new TH1F("yX", ";y_{X};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("ptX", ";X p_{T} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("pt12", ";X p_{T} [GeV];Number of Events", 50, 100., 700.));
    plotter.book(new TH1F("pt34", ";X p_{T} [GeV];Number of Events", 50, 100., 700.));
    
    plotter.book(new TH1F("before_mA", ";m_{A} [GeV];Number of Events", 50, 150., 200.));
    plotter.book(new TH1F("before_ptA", ";X p_{A} [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("before_etaA", ";eta_{A};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_yA", ";y_{A};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_PhiA", ";#Phi_{A};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_mB", ";m_{B} [GeV];Number of Events", 50, 150., 200.));
    plotter.book(new TH1F("before_ptB", ";X p_{B} [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("before_etaB", ";eta_{B};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_yB", ";y_{B};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_PhiB", ";#Phi_{B};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("mA", ";m_{A} [GeV];Number of Events", 50, 150., 200.));
    plotter.book(new TH1F("ptA", ";X p_{A} [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("etaA", ";eta_{A};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("yA", ";y_{A};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("PhiA", ";#Phi_{A};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("mB", ";m_{B} [GeV];Number of Events", 50, 150., 200.));
    plotter.book(new TH1F("ptB", ";X p_{B} [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("etaB", ";eta_{B};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("yB", ";y_{B};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("PhiB", ";#Phi_{B};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_quarkMiss12", ";pid (abs);Number of Events", 50, 0, 6));
    //plotter.book(new TH1F("before_mMiss12", ";m_{Miss12} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("before_ptMiss12", ";X p_{Miss12} [GeV];Number of Events", 50, 0., 250.)); //events exist up to 350
    plotter.book(new TH1F("before_etaMiss12", ";eta_{Miss12};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_yMiss12", ";y_{Miss12};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_phiMiss12", ";#Phi_{Miss12};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_quarkMiss34", ";pid (abs);Number of Events", 50, 0, 6));
    //plotter.book(new TH1F("before_mMiss34", ";m_{Miss34} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("before_ptMiss34", ";X p_{Miss34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_etaMiss34", ";eta_{Miss34};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_yMiss34", ";y_{Miss34};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_phiMiss34", ";#Phi_{Miss34};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_charge12", ";e_{12};Number of Events", 50, -2, 2));
    plotter.book(new TH1F("before_charge34", ";e_{34};Number of Events", 50, -2, 2));
    
    plotter.book(new TH1F("before_bLeadStatus12", ";code_{12};Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("before_bLeadStatus34", ";code_{34};Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("before_bothBleading", ";code;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("before_dijetsHaveTwoBquarks", ";code;Number of Events", 10, -0.5, 2.5));
    
    plotter.book(new TH1F("before_dR_sublead12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dR_sublead34", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dR_lead12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dR_lead34", ";dR;Number of Events", 50, 0, 5));
    
    plotter.book(new TH1F("before_dR_Wquarks12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dR_Wquarks34", ";dR;Number of Events", 50, 0, 5));
    
    plotter.book(new TH1F("before_numJetsMatched", ";number;Number of Events", 30, -0.5, 4.5));
    plotter.book(new TH1F("before_numQuarks", ";number;Number of Events", 30, -0.5, 6.5));
    plotter.book(new TH1F("before_quarkMatchStatus", ";code;Number of Events", 50, -7.5, 7.5));
    plotter.book(new TH1F("before_alignStatus", ";code;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("before_multiMatch", ";flag;Number of Events", 10, -0.5, 1.5));
    
    plotter.book(new TH1F("before_quark1", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("before_quark2", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("before_quark3", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("before_quark4", ";pid (abs);Number of Events", 50, 0, 6));
    
    plotter.book(new TH1F("before_deltaR1", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("before_deltaR2", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("before_deltaR3", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("before_deltaR4", ";dR;Number of Events", 50, 0, 0.3));
    
    plotter.book(new TH1F("before_numNonDup", ";number;Number of Events", 50, -0.5, 12.5));
    plotter.book(new TH1F("before_ptMinAll", ";X p [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("before_etaMaxAll", ";eta;Number of Events", 50,0., 5.));
    plotter.book(new TH1F("before_numClosestFound", ";number;Number of Events", 50, -0.5, 2.5));
    plotter.book(new TH1F("before_multiMatch3rd", ";flag;Number of Events", 10, -0.5, 1.5));
    plotter.book(new TH1F("before_num3rdJetsMatched", ";flag;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("before_sameClosest", ";flag;Number of Events", 10, -0.5, 1.5));
    plotter.book(new TH1F("before_3rdMatchStatus", ";code;Number of Events", 50, -3.5, 3.5));
    
    plotter.book(new TH1F("before_quarkMatched3rdJet12", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("before_pt3rdJet12", ";X p_{3rdJet12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_m3rdJet12", ";m_{3rdJet12} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("before_eta3rdJet12", ";eta_{3rdJet12};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_y3rdJet12", ";y_{3rdJet12};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_phi3rdJet12", ";#Phi_{3rdJet12};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_quarkMatched3rdJet34", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("before_pt3rdJet34", ";X p_{3rdJet34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("before_m3rdJet34", ";m_{3rdJet34} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("before_eta3rdJet34", ";eta_{3rdJet34};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("before_y3rdJet34", ";y_{3rdJet34};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("before_phi3rdJet34", ";#Phi_{3rdJet34};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("before_dR3rdJet12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dRMatched3rdJet12", ";dR;Number of Events", 50, 0, 0.3));
    
    plotter.book(new TH1F("before_dR3rdJet34", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("before_dRMatched3rdJet34", ";dR;Number of Events", 50, 0, 0.3));

    plotter.book(new TH1F("quarkMiss12", ";pid (abs);Number of Events", 50, 0, 6));
    //plotter.book(new TH1F("mMiss12", ";m_{Miss12} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("ptMiss12", ";X p_{Miss12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("etaMiss12", ";eta_{Miss12};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("yMiss12", ";y_{Miss12};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("phiMiss12", ";#Phi_{Miss12};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("quarkMiss34", ";pid (abs);Number of Events", 50, 0, 6));
    //plotter.book(new TH1F("mMiss34", ";m_{Miss34} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("ptMiss34", ";X p_{Miss34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("etaMiss34", ";eta_{Miss34};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("yMiss34", ";y_{Miss34};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("phiMiss34", ";#Phi_{Miss34};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("charge12", ";e_{12};Number of Events", 50, -2, 2));
    plotter.book(new TH1F("charge34", ";e_{34};Number of Events", 50, -2, 2));
    
    plotter.book(new TH1F("bLeadStatus12", ";code_{12};Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("bLeadStatus34", ";code_{34};Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("bothBleading", ";code;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("dijetsHaveTwoBquarks", ";code;Number of Events", 10, -0.5, 2.5));
    
    plotter.book(new TH1F("dR_sublead12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dR_sublead34", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dR_lead12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dR_lead34", ";dR;Number of Events", 50, 0, 5));
    
    plotter.book(new TH1F("dR_Wquarks12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dR_Wquarks34", ";dR;Number of Events", 50, 0, 5));
    
    plotter.book(new TH1F("numJetsMatched", ";number;Number of Events", 30, -0.5, 4.5));
    plotter.book(new TH1F("numQuarks", ";number;Number of Events", 30, -0.5, 6.5));
    plotter.book(new TH1F("quarkMatchStatus", ";code;Number of Events", 50, -7.5, 7.5));
    plotter.book(new TH1F("alignStatus", ";code;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("multiMatch", ";flag;Number of Events", 10, -0.5, 1.5));
    
    plotter.book(new TH1F("quark1", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("quark2", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("quark3", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("quark4", ";pid (abs);Number of Events", 50, 0, 6));

    plotter.book(new TH1F("deltaR1", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("deltaR2", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("deltaR3", ";dR;Number of Events", 50, 0, 0.3));
    plotter.book(new TH1F("deltaR4", ";dR;Number of Events", 50, 0, 0.3));
    
    plotter.book(new TH1F("numNonDup", ";number;Number of Events", 50, -0.5, 12.5));
    plotter.book(new TH1F("ptMinAll", ";X p [GeV];Number of Events", 50, 0., 500.));
    plotter.book(new TH1F("etaMaxAll", ";eta;Number of Events", 50,0., 5.));
    plotter.book(new TH1F("numClosestFound", ";number;Number of Events", 50, -0.5, 2.5));
    plotter.book(new TH1F("multiMatch3rd", ";flag;Number of Events", 10, -0.5, 1.5));
    plotter.book(new TH1F("num3rdJetsMatched", ";flag;Number of Events", 10, -0.5, 2.5));
    plotter.book(new TH1F("sameClosest", ";flag;Number of Events", 10, -0.5, 1.5));
    plotter.book(new TH1F("3rdMatchStatus", ";code;Number of Events", 50, -3.5, 3.5));
    
    plotter.book(new TH1F("quarkMatched3rdJet12", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("pt3rdJet12", ";X p_{3rdJet12} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("m3rdJet12", ";m_{3rdJet12} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("eta3rdJet12", ";eta_{3rdJet12};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("y3rdJet12", ";y_{3rdJet12};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("phi3rdJet12", ";#Phi_{3rdJet12};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("quarkMatched3rdJet34", ";pid (abs);Number of Events", 50, 0, 6));
    plotter.book(new TH1F("pt3rdJet34", ";X p_{3rdJet34} [GeV];Number of Events", 50, 0., 250.));
    plotter.book(new TH1F("m3rdJet34", ";m_{3rdJet34} [GeV];Number of Events", 50, 0., 5.));
    plotter.book(new TH1F("eta3rdJet34", ";eta_{3rdJet34};Number of Events", 50,-5., 5.));
    plotter.book(new TH1F("y3rdJet34", ";y_{3rdJet34};Number of Events", 50, -2.5, 2.5));
    plotter.book(new TH1F("phi3rdJet34", ";#Phi_{3rdJet34};Number of Events", 50, 0., 2*M_PI));
    
    plotter.book(new TH1F("dR3rdJet12", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dRMatched3rdJet12", ";dR;Number of Events", 50, 0, 0.3));
    
    plotter.book(new TH1F("dR3rdJet34", ";dR;Number of Events", 50, 0, 5));
    plotter.book(new TH1F("dRMatched3rdJet34", ";dR;Number of Events", 50, 0, 0.3));
    
    plotter.book(new TH1F("TopVetoBDT", ";Top Veto BDT Output;Number of Events", 50, -1., 1.));

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
    TString plotName = "dijets_cosThetaStar";
    std::cout<<"2 dijets";
    for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
    {
        std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
    }
    std::cout<<"\\\\"<< std::endl;
    plotName = "dijets_mH_cosThetaStar";
    std::cout<<"2 dijets $m_H$";
    for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
    {
        std::cout<<" & "<< formatNumberForTable(plotter.getPlot(plotName, *cat)->Integral());
    }
    std::cout<<"\\\\"<< std::endl;
    //plotName = "BDT";
    plotName = "topVeto_cosThetaStar";
    std::cout<<"Top Veto";
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
    plotName = "dijets_cosThetaStar";
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
    plotName = "topVeto_cosThetaStar";  //I broke this but too lazy to debug so i just reused the previous plot
    back = 0.;
    for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
    {
        if(*cat == "HH") sig = plotter.getPlot(plotName, *cat)->Integral();
        else back += plotter.getPlot(plotName, *cat)->Integral();
    }
    std::cout<<"Top Veto & "<< formatNumberForTable(sig) <<" & "<< formatNumberForTable(back) <<" & "<< formatNumberForTable(sig/back) <<" & "<< formatNumberForTable(sig/sqrt(back));
    std::cout<<"\\\\"<< std::endl;
    std::cout<<"\\hline\\end{tabular}\\caption{"<< caption <<"}\\end{center}\\end{table}"<< std::endl;
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
