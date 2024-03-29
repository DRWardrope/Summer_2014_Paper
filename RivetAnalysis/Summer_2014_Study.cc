// -*- C++ -*-
#include <algorithm>
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "FinalStateWithGhosts.hh"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TagJet.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


	class Summer_2014_Study : public Analysis 
	{
		public:

			/// Constructor
			Summer_2014_Study()
				: Analysis("Summer_2014_Study")
			{    
			}


		public:

			/// Book histograms and initialise projections before the run
			void init()
			{
				//_crossSection = handler().crossSection();
				//Create final State projection with added ghost particles.
				FinalStateWithGhosts fs(-5.0, 5.0);
				fs.ghostIdPair(PID::BQUARK);
				fs.ghostIdPair(PID::CQUARK);
				fs.ghostIdPair(PID::TAU);
				addProjection(fs, "FS");

				FastJets akt4(fs, FastJets::ANTIKT, 0.4);
				addProjection(akt4, "AntiKt4");

				bookTree();

			}


			/// Perform the per-event analysis
			void analyze(const Event& event) 
			{
				const double weight = event.weight();
				bookAndFill1D("CutFlow", "Yields vs Selection Step", 0, "Selection Step", 5, -0.5, 4.5, weight);

				const Jets &akt4 = applyProjection<JetAlg>(event, "AntiKt4").jetsByPt(40*GeV);
				std::vector<TagJet> bjets;
				bjets.reserve(4);
				foreach(const Jet &j, akt4) 
				{
					if (fabs(j.momentum().rapidity()) < 2.5)
					{
						if(j.containsBottom())
						{
							bjets.push_back(TagJet(j, 0.7, 5));
						}else if(j.containsCharm()) 
						{
							bjets.push_back(TagJet(j, 0.2, 4));
						}else if(j.containsParticleId(15) || j.containsParticleId(-15))
						{
							bjets.push_back(TagJet(j, 0.2, 15));
						}else 
						{
							bjets.push_back(TagJet(j, 0.01, 0));
						}
					}
				}
				//Demand at least 4 jets with pt > 40 GeV and |eta| < 2.5. Despite name, not yet b-tagged!
				if(bjets.size() < 4) return;
				bookAndFill1D("CutFlow", "Yields vs Selection Step", 1, "Selection Step", 5, -0.5, 4.5, weight);

				if (dijetSelection(bjets)) 
				{
					FourMomentum dijet1 = bjets[0].momentum()+bjets[1].momentum();
					FourMomentum dijet2 = bjets[2].momentum()+bjets[3].momentum();
					const float btaggedWeight = weight*bjets[0].tagEff()*bjets[1].tagEff()*bjets[2].tagEff()*bjets[3].tagEff();

					fillKinematics("dijet", "Dijet", dijet1, btaggedWeight);
					fillKinematics("dijet", "Dijet", dijet2, btaggedWeight);
					for (unsigned int i=0; i<4; ++i) {
						fillKinematics("jet", "Jet", bjets[i].momentum(), btaggedWeight);
					}
					fillXttVariables(akt4,"jet", "Jet", dijet1, bjets[0], bjets[1], dijet2, bjets[2], bjets[3], btaggedWeight); 
					fillMelaAngles(event, "jet", "Jet", dijet1, bjets[0], bjets[1], dijet2, bjets[2], bjets[3], btaggedWeight);
					//dijetSelection only satisfied if ≥2 jets are found. Use weight, since x=2 represents 2 dijets w/o b-tagging req.
					bookAndFill1D("CutFlow", "Yields vs Selection Step", 2, "Selection Step", 5, -0.5, 4.5, weight);
					//Use btaggedWeight to replicate effect of b-tagging.
					bookAndFill1D("CutFlow", "Yields vs Selection Step", 3, "Selection Step", 5, -0.5, 4.5, btaggedWeight);
				}

			}


			/// Normalise histograms etc., after the run
			void finalize() 
			{
				int nEvent_input = -99;
				double sumW_input = -99.;
				std::map<std::string, boost::shared_ptr<YODA::Histo1D> >::const_iterator hIt = _histograms_1d.find("CutFlow");
				if(hIt == _histograms_1d.end()) std::cout <<"finalize: Could not obtain CutFlow YODA::Histo1D!"<< std::endl;
				else{
					const YODA::HistoBin1D& bin0 = hIt->second->bin(hIt->second->binIndexAt(0));
					sumW_input = bin0.area();
					nEvent_input = bin0.numEntries();
				}
				double xSection = crossSection();
				double intLumi = nEvent_input/xSection;

				TTree* tRunInfo = new TTree("RunInfo", "Aggregated information for each analysis job");
				tRunInfo->Branch("nEvent_input", &nEvent_input, "nEvent_input/I");
				tRunInfo->Branch("sumW_input", &sumW_input, "sumW_input/D");
				tRunInfo->Branch("crossSection", &xSection, "crossSection/D");
				tRunInfo->Branch("intLumi", &intLumi, "intLumi/D");
				tRunInfo->Fill();

				//std::cout<<"In the analysis, we have found n_events = "<< 

				_fOut->Write();
				_fOut->Close();
				/// @todo Normalise, scale and otherwise manipulate histograms here

			/*	map<string, boost::shared_ptr<YODA::Histo1D> >::iterator it1;
				for (it1=_histograms_1d.begin(); it1!=_histograms_1d.end(); ++it1) {
					normalize(it1->second);
				}*/

				/*map<string, boost::shared_ptr<YODA::Histo2D> >::iterator it2;
				  for (it2=_histograms_2d.begin(); it2!=_histograms_2d.end(); ++it2) {
				// Rivet crashes if you don't normalize your 2D histograms (???)
				normalize(it2->second);
				}*/

				// Create efficiency histograms
				finalizeEfficiency("akt4_H_pt_", "separate", "merged", "Higgs $p_T$ [GeV]");
				finalizeEfficiency("akt4_bb_dr_", "separate", "merged", "$\\Delta R(b,b)$");
				finalizeEfficiency("ca12_H_pt_", "merged", "separate", "Higgs $p_T$ [GeV]");
				finalizeEfficiency("ca12_bb_dr_", "merged", "separate", "$\\Delta R(b,b)$");
				/// @todo Normalise, scale and otherwise manipulate histograms here

				// scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
				// normalize(_h_YYYY); // normalize to unity

			}
			//This sets the structure of the output TTree and associates the TBranches with variables.
			void bookTree()
			{
				_fOut = TFile::Open("MVAResults.root", "RECREATE");
				/*if((crossSection() - 146.03) < 0.0001) _tree = new TTree("bbbb", "TMVA Input Tree");
				else if((crossSection() - 317.78) < 0.0001) _tree = new TTree("bbcc", "TMVA Input Tree");
				else if((crossSection() - 212.07) < 0.0001) _tree = new TTree("tt", "TMVA Input Tree");
				else{	*/
				std::cout<<"bookTree: couldn't figure out a TTree name based on the magic cross-section numbers."<<std::endl;
				_tree = new TTree("TMVAInput", "TMVA Input Tree");
					//_tree = new TTree(handler().runName().c_str(), "TMVA Input Tree");
				//}				
				_tree->Branch("weight", &_weight, "weight/D");
				_tree->Branch("btagWeight", &_btagWeight, "btagWeight/D");
				_tree->Branch("genWeight", &_genWeight, "genWeight/D");

				_tree->Branch("mX", &_mX, "mX/D");
				_tree->Branch("yX", &_yX, "yX/D");
				_tree->Branch("ptX", &_ptX, "ptX/D");
				_tree->Branch("etaX", &_etaX, "etaX/D");
				_tree->Branch("phiX", &_phiX, "phiX/D");

				_tree->Branch("pt12", &_pt12, "pt12/D");
				_tree->Branch("eta12", &_eta12, "eta12/D");
				_tree->Branch("phi12", &_phi12, "phi12/D");
				_tree->Branch("m12", &_m12, "m12/D");
				_tree->Branch("pt34", &_pt34, "pt34/D");
				_tree->Branch("eta34", &_eta34, "eta34/D");
				_tree->Branch("phi34", &_phi34, "phi34/D");
				_tree->Branch("m34", &_m34, "m34/D");
	
				_tree->Branch("mW12", &_mW12, "mW12/D");
				_tree->Branch("mt12", &_mt12, "mt12/D");
				_tree->Branch("dRW12", &_dRW12, "dRW12/D");
				_tree->Branch("mW34", &_mW34, "mW34/D");
				_tree->Branch("mt34", &_mt34, "mt34/D");
				_tree->Branch("dRW34", &_dRW34, "dRW34/D");

				_tree->Branch("pt1", &_pt1, "pt1/D");
				_tree->Branch("eta1", &_eta1, "eta1/D");
				_tree->Branch("phi1", &_phi1, "phi1/D");
				_tree->Branch("m1", &_m1, "m1/D");
				_tree->Branch("flav1", &_flav1, "flav1/I");
				_tree->Branch("pt2", &_pt2, "pt2/D");
				_tree->Branch("eta2", &_eta2, "eta2/D");
				_tree->Branch("phi2", &_phi2, "phi2/D");
				_tree->Branch("m2", &_m2, "m2/D");
				_tree->Branch("flav2", &_flav2, "flav2/I");
				_tree->Branch("pt3", &_pt3, "pt3/D");
				_tree->Branch("eta3", &_eta3, "eta3/D");
				_tree->Branch("phi3", &_phi3, "phi3/D");
				_tree->Branch("m3", &_m3, "m3/D");
				_tree->Branch("flav3", &_flav3, "flav3/I");
				_tree->Branch("pt4", &_pt4, "pt4/D");
				_tree->Branch("eta4", &_eta4, "eta4/D");
				_tree->Branch("phi4", &_phi4, "phi4/D");
				_tree->Branch("m4", &_m4, "m4/D");
				_tree->Branch("flav4", &_flav4, "flav4/I");

				_tree->Branch("cosThetaStar", &_cosThetaStar, "cosThetaStar/D");
				_tree->Branch("Phi", &_Phi, "Phi/D");
				_tree->Branch("Phi1", &_Phi1, "Phi1/D");
				_tree->Branch("cosTheta1", &_cosTheta1, "cosTheta1/D");
				_tree->Branch("cosTheta2", &_cosTheta2, "cosTheta2/D");

			}
			//fillMelaAngles fills the output TTree with angles calculated by calculateAngles().
			void fillMelaAngles(const Event& event, std::string label, std::string title, const FourMomentum &H1, const TagJet &b11, 
					const TagJet &b12, const FourMomentum &H2, const TagJet &b21, const TagJet &b22, double weight) 
			{

				double cos_theta_star, phi, phi1, cos_theta1, cos_theta2;
				calculateAngles(b11.momentum(), b12.momentum(), b21.momentum(), b22.momentum(), cos_theta_star, phi, phi1, cos_theta1, cos_theta2);
				FourMomentum X = H1 + H2;
				if(label == "jet")
				{
					_genWeight = event.weight();
					_weight = weight;
					_btagWeight = b11.tagEff()*b12.tagEff()*b21.tagEff()*b22.tagEff();
					_cosThetaStar = cos_theta_star;
					_cosTheta1 = cos_theta1;
					_cosTheta2 = cos_theta2;
					_Phi = phi;
					_Phi1 = phi1;
					//Properties of Higgs bosons
					_pt12 = H1.pt();
					_eta12 = H1.eta();
					_phi12 = H1.phi();
					_m12 = H1.mass();
					_pt34 = H2.pt();
					_eta34 = H2.eta();
					_phi34 = H2.phi();
					_m34 = H2.mass();
					//X properties
					_mX = X.mass();
					_yX = X.rapidity();
					_etaX = X.eta();
					_phiX = X.phi();
					_ptX = X.pt();
					//Individual jet properties
					_pt1 = b11.pt();
					_eta1 = b11.eta();
					_phi1 = b11.phi();
					_m1 = b11.mass();
					_flav1 = b11.flavour();
					_pt2 = b12.pt();
					_eta2 = b12.eta();
					_phi2 = b12.phi();
					_m2 = b12.mass();
					_flav2 = b12.flavour();
					_pt3 = b21.pt();
					_eta3 = b21.eta();
					_phi3 = b21.phi();
					_m3 = b21.mass();
					_flav3 = b21.flavour();
					_pt4 = b22.pt();
					_eta4 = b22.eta();
					_phi4 = b22.phi();
					_m4 = b22.mass();
					_flav4 = b22.flavour();
	
					//std::cout<<"Filling TTree with MELA information."<< std::endl;
					_tree->Fill();
				}	
				bookAndFill1D(label+"_mX", title+" $m_{X}$",   X.mass(),   "$m_{X} [GeV]$",   50, 0.,   1000.,   weight);
				bookAndFill1D(label+"_ptX", title+" $pt_{X}$",   X.pt(),   "$X p_{T} [GeV]$",   50, 0.,   500.,   weight);
				bookAndFill1D(label+"_yX", title+" $y_{X}$",   X.rapidity(),   "$y_X$",   50, -2.,   2.,   weight);
				bookAndFill1D(label+"_m12", title+" $m_{12}$",   H1.mass(),   "$m_{12} [GeV]$",   50, 0.,   200.,   weight);
				bookAndFill1D(label+"_m34", title+" $m_{34}$",   H2.mass(),   "$m_{34} [GeV]$",   50, 0.,   200.,   weight);
				bookAndFill1D(label+"_cos_theta_star", title+" $\\cos(\\theta^*)$",   cos_theta_star,   "$\\cos(\\theta^*)$",   50, -1.,   1.,   weight);
				bookAndFill1D(label+"_phi",            title+" $\\Phi$",              phi,              "$\\Phi$",              50, -M_PI, M_PI, weight);
				bookAndFill1D(label+"_phi1",           title+" $\\Phi_1$",            phi,              "$\\Phi_1$",            50, -M_PI, M_PI, weight);
				bookAndFill1D(label+"_cos_theta1",     title+" $\\cos(\\theta_1)$",   cos_theta1,       "$\\cos(\\theta_1)$",   50, -1.,   1.,   weight);
				bookAndFill1D(label+"_cos_theta2",     title+" $\\cos(\\theta_2)$",   cos_theta2,       "$\\cos(\\theta_2)$",   50, -1.,   1.,   weight);
				bookAndFill1D(label+"_abs_cos_theta1", title+" $|\\cos(\\theta_1)|$", fabs(cos_theta1), "$|\\cos(\\theta_1)|$", 50,  0.,   1.,   weight);
				bookAndFill1D(label+"_abs_cos_theta2", title+" $|\\cos(\\theta_2)|$", fabs(cos_theta2), "$|\\cos(\\theta_2)|$", 50,  0.,   1.,   weight);
			}

			//fillXttVariables calculates variables that are useful for rejecting top-quark backgrounds. 
			void fillXttVariables(const Jets& jets, std::string label, std::string title, const FourMomentum &H1, const TagJet &b11, 
					const TagJet &b12, const FourMomentum &H2, const TagJet &b21, const TagJet &b22, double weight) 
			{
				_mW12 = -99, _mW34 = -99, _mt12 = -99, _mt34 = -99, _dRW12 = -99, _dRW34 = -99;
				Jets noDups;
				foreach(const Jet& j, jets)
				{
					if(j.eta() == b11.eta() && j.phi() == b11.phi()) continue;
					if(j.eta() == b12.eta() && j.phi() == b12.phi()) continue;
					if(j.eta() == b21.eta() && j.phi() == b21.phi()) continue;
					if(j.eta() == b22.eta() && j.phi() == b22.phi()) continue;
					noDups.push_back(j);
				}
				assert(noDups.size() == jets.size() - 4);
				const TagJet &leastTagged12 = b11.tagEff() >= b12.tagEff() ? b12 : b11; 
				const TagJet &leastTagged34 = b21.tagEff() >= b22.tagEff() ? b22 : b21; 
				//float mindR12 = 1.5, mindR34 = 1.5;
				float mindR12 = 2.0, mindR34 = 2.0;
				const Jet *closestJet12=NULL, *closestJet34=NULL;
				foreach(const Jet& j, noDups)
				{
					float dR12 = deltaR(leastTagged12, j);			
					if(dR12 < mindR12 && dR12 > 0.1) 
					{
						mindR12 = dR12;
						closestJet12 = &j;
					}
					float dR34 = deltaR(leastTagged34, j);			
					if(dR34 < mindR34 && dR34 > 0.1)
					{
						mindR34 = dR34;
						closestJet34 = &j;
					}
				}
				if(closestJet12)
				{
					_mW12 = (closestJet12->momentum() + leastTagged12.momentum()).mass();
					_mt12 = (closestJet12->momentum() + b11.momentum() + b12.momentum()).mass();	
					_dRW12 = mindR12;
				}
				if(closestJet34)
				{
					_mW34 = (closestJet34->momentum() + leastTagged34.momentum()).mass();
					_mt34 = (closestJet34->momentum() + b21.momentum() + b22.momentum()).mass();	
					_dRW34 = mindR34;
				}
			}

			/// Fill a histogram with the unsigned PDG ID of daughter particles for given parent particle
			void fillDaughters(std::string label, std::string title, const HepMC::GenParticle *parent, double weight) {

				const HepMC::GenVertex *vtx = parent->end_vertex();
				if (!vtx) return;

				HepMC::GenVertex::particles_out_const_iterator it;
				for (it=vtx->particles_out_const_begin(); it!=vtx->particles_out_const_end(); ++it) {
					const HepMC::GenParticle *daughter = *it;
					bookAndFill1D(label+"_daughters", title, abs(daughter->pdg_id()), "Unsigned PDG ID", 40,  0, 40, weight);
				}
			}

			/// Fill histograms with angular variables dphi, deta, dr for the two four vectors
			void fillAngles(std::string label, std::string title, const FourMomentum &momentum1, const FourMomentum &momentum2, double weight) {

				double dphi = deltaPhi(momentum1.phi(), momentum2.phi());
				double deta = deltaEta(momentum1.eta(), momentum2.eta());
				double dr   = hypot(dphi, deta);

				bookAndFill1D(label+"_dphi", title, dphi, "$\\Delta\\phi$", 100, 0., M_PI, weight);
				bookAndFill1D(label+"_deta", title, dphi, "$\\Delta\\eta$", 100, 0., 4.,   weight);
				bookAndFill1D(label+"_dr",   title, dr,   "$\\Delta R$",    100, 0., 6.,   weight);
			}

			/// Fill some kinematic variable histograms for the given four vector
			void fillKinematics(std::string label, std::string title, const FourMomentum &momentum, double weight) {

				bookAndFill1D(label+"_pt",   title, momentum.perp()*GeV, "$p_T$ [GeV]", 200,  0., 1e3,     weight);
				bookAndFill1D(label+"_eta",  title, momentum.eta(),      "$\\eta$",     100, -5., 5.,      weight);
				bookAndFill1D(label+"_phi",  title, momentum.phi(),      "$\\phi$",     100,  0., 2.*M_PI, weight);

				double max = label.compare("BB") ? 1e3 : 3e3; // 1TeV normally, but for the BB pair 3TeV
				bookAndFill1D(label+"_mass", title, momentum.mass()*GeV, "$m$ [GeV]",   200,  0., max,     weight);
			}

			void fillImpactParameters(std::string label, std::string title, const GenVertex *vertex, double weight) {

				if (!vertex) return;
				if (vertex->point3d().perp() == 0. && vertex->point3d().z() == 0.) return;


				bookAndFill1D(label+"_d0",   title, vertex->point3d().perp(), "$|d_0|$ [mm]", 100, 0.,    100,   weight);
				bookAndFill1D(label+"_z0",   title, vertex->point3d().z(),    "$z_0$ [mm]",   100, -100., 100.,  weight);
				bookAndFill1D(label+"_phi0", title, vertex->point3d().phi(),  "$\\phi_0$",    100,  -M_PI, M_PI, weight);
			}

			const HepMC::GenParticle *findLastInChain(const HepMC::GenParticle *particle) {
				HepMC::GenVertex *decay = particle->end_vertex();
				if (decay) {
					HepMC::GenVertex::particles_out_const_iterator it;
					for (it=decay->particles_out_const_begin(); it!=decay->particles_out_const_end(); ++it) {
						if ((*it)->pdg_id() == particle->pdg_id() && (*it)->barcode() > particle->barcode()) {
							return findLastInChain(*it);
						}
					}
				}
				return particle;
			}

			void fillReconstructionEfficiency(std::string label, std::string title, const Particle &parent, 
					const Particle &p1, const Particle &p2, const Jets &jets, double weight) {

				const Jet *j1 = 0, *j2 = 0;
				foreach(const Jet &j, jets) {
					if (j.containsParticle(p1)) {
						j1 = &j;
					} 
					if (j.containsParticle(p2)) {
						j2 = &j;
					}
				}

				std::string result;
				if (j1 && j2) {
					result = j1 == j2 ? "merged" : "separate";
				} else {
					result = "missed";
				}

				bookAndFill1D(label+"_H_pt_"+result,  title, parent.momentum().pT()/GeV,           "Higgs $p_T$ [GeV]", 60, 150., 750., weight);
				bookAndFill1D(label+"_bb_dr_"+result, title, deltaR(p1.momentum(), p2.momentum()), "$\\Delta R(b,b)$",  50, 0.,   2.0,  weight);
			}

			void finalizeEfficiency(std::string label, std::string pass, std::string fail, std::string axis="") {
/*				boost::shared_ptr<YODA::Histo1D> hist_pass = _histograms_1d[label+pass];
				boost::shared_ptr<YODA::Histo1D> hist_fail = _histograms_1d[label+fail];

				if (!hist_pass) return;

				for (unsigned int i=0; i<hist_pass->bins().size(); ++i) {
					double entries_pass = hist_pass->bin(i).sumW();
					double entries_total = entries_pass + (hist_fail ? hist_fail->bin(i).sumW() : 0.);

//					if (entries_total > 0 && entries_pass >= 0) {
//					  bookAndFill1D(label+"efficiency", hist_pass->title(), (hist_pass->bin(i).lowEdge()+hist_pass->bin(i).highEdge())/2., 
//					  axis, std::string("BLARG"), hist_pass->lowEdge(), hist_pass->highEdge(), 
//
//					//            axis, hist_pass->bins()->title(), hist_pass->lowEdge(), hist_pass->highEdge(), 
//					entries_pass/entries_total * hist_pass->bin(i).width());
//					}
				}*/
			}

			/// Histogram booking and filling in one step, histograms are automatically created at the first fill
			void bookAndFill1D(std::string name, std::string title, 
					double x, std::string xaxis, unsigned int nbinsx, double xlow, double xup, 
					double weight) 
			{
				if(isnan(x)) 
				{
					std::cout <<"bookAndFill1D: x = "<< x <<", returning."<< std::endl;
					return;
				}
				//std::cout<<"bookAndFill1D("<< name <<", "<< title <<", "<< x <<")"<< std::endl;
				boost::shared_ptr<YODA::Histo1D> hist = _histograms_1d[name];
				if (!hist) {
					std::cout<<"bookAndFill1D: "<< name <<" undefined, book a new one with [" << xlow <<", "<< xup <<"]"<< std::endl;
					hist = bookHisto1D(name, nbinsx, xlow, xup, title, xaxis);
					_histograms_1d[name] = hist;
				}
				hist->fill(x, weight);
			}

			void bookAndFill2D(std::string name, std::string title, 
					double x, std::string xaxis, unsigned int nbinsx, double xlow, double xup, 
					double y, std::string yaxis, unsigned int nbinsy, double ylow, double yup, 
					double weight) {

				boost::shared_ptr<YODA::Histo2D> hist = _histograms_2d[name];
				if (!hist) {
					//hist = bookHisto2D(name, nbinsx, xlow, xup, nbinsy, ylow, yup, title, xaxis, yaxis);
					// _histograms_2d[name] = hist;
				}
				//hist->fill(x, y, weight);
			}

			//tagGreater is used when sorting vector<TagJet> by b-tagging weight
        	static bool tagGreater(const TagJet& a, const TagJet& b){ return (a.tagEff() > b.tagEff()); }
			//dijetSelection tries to find Higgs boson candidates from an input of vector<TagJet>.
			//It re-orders the jets in this vector, such that bjets[0] and [1] correspond to Higgs 1 and [2] and [3] to Higgs 2
			bool dijetSelection(std::vector<TagJet> &bjets) 
			{
				std::vector<TagJet> input = bjets;	
				std::stable_sort(input.begin(), input.end(), tagGreater);
				/*for(int j = 0; j < input.size(); ++j)
				{
					std::cout<<"dijetSelection: "<< j <<", tagEff = "<< input[j].tagEff() <<", pt = "<< input[j].pt() << std::endl;
				}*/
				bjets.clear();

				/*std::cout<<"Input jet collection:"<<std::endl;
				for (unsigned int i=0; i<input.size(); ++i) 
				{
					std::cout<< i <<": "<< input[i].pt() <<" GeV"<< std::endl;
				}
				//std::cout <<"------------------------------------------------"<< std::endl;*/
				for (unsigned int i=0; i<input.size(); ++i) 
				{
					//std::cout<<"jet "<< i <<" pt = "<< input[i].pt() << std::endl;
					for (unsigned int j=i+1; j<input.size(); ++j) 
					{
						if (deltaR(input[i], input[j]) > 1.5) continue;
						FourMomentum dijet = input[i].momentum() + input[j].momentum();

						if (dijet.pT() < 150*GeV) continue;
						//if (fabs(dijet.mass()/GeV - 115) > 25) continue;
					//	std::cout<<"Dijet created from inputs "<<i <<" and "<< j <<". pt1 = "<< input[i].pt() <<", pt2 = "<< input[j].pt() <<" GeV"<< std::endl;
						bjets.push_back(input[i]);
						bjets.push_back(input[j]);

							
						input.erase(input.begin() + j--);
						input.erase(input.begin() + i--);//erases element at begin()+i, then decrements i by one, so we don't miss out an element because of the vector shortening.
						/*std::cout<<"Input jet collection now contains:"<<std::endl;
						for (unsigned int k=0; k<input.size(); ++k) 
						{
							std::cout<< k <<": "<< input[k].pt() <<" GeV"<< std::endl;
						}
						std::cout <<"..............................."<< std::endl;*/
						break;
					}
				}

				return bjets.size() >= 4;
			}

			//calculateAngles works out the decay angles of the Higgs bosons.
			void calculateAngles(const FourMomentum &q11_lab, const FourMomentum &q12_lab, 
					const FourMomentum &q21_lab, const FourMomentum &q22_lab,
					double &cos_theta_star, double &phi, double &phi1, 
					double &cos_theta1, double &cos_theta2) {

				FourMomentum q1_lab = q11_lab + q12_lab;
				FourMomentum q2_lab = q21_lab + q22_lab;
				FourMomentum X_lab = q1_lab + q2_lab;

				LorentzTransform to_X_CoM(-X_lab.boostVector());
				FourMomentum q1  = to_X_CoM.transform(q1_lab);
				FourMomentum q11 = to_X_CoM.transform(q11_lab);
				FourMomentum q12 = to_X_CoM.transform(q12_lab);
				FourMomentum q2  = to_X_CoM.transform(q2_lab);
				FourMomentum q21 = to_X_CoM.transform(q21_lab);
				FourMomentum q22 = to_X_CoM.transform(q22_lab);

				//PRINT(q1_lab);
				//PRINT(q2_lab);
				//PRINT(q1_lab+q2_lab);
				//PRINT(q1);
				//PRINT(q2);
				//PRINT(q1+q2);

				assert(fuzzyEquals((q1+q2).pT(), 0., 1e-4)); 
				assert(fuzzyEquals(q1, q11+q12, 1e-4));
				assert(fuzzyEquals(q2, q21+q22, 1e-4));

				Vector3 nZ(0., 0., 1.);
				Vector3 n1  = q11.p3().cross(q12.p3()).unit();
				Vector3 n2  = q21.p3().cross(q22.p3()).unit();
				Vector3 nSC = nZ.cross(q1.p3()).unit();

				//PRINT(n1);
				//PRINT(n2);
				//PRINT(nSC);

				cos_theta_star = cos(q1.p3().theta()-nZ.theta());
				phi  = sign(q1.p3().dot(n1.cross(n2))) * acos(-n1.dot(n2)); 
				phi1 = sign(q1.p3().dot(n1.cross(nSC))) * acos(n1.dot(nSC));

				LorentzTransform to_H1_CoM(-q1_lab.boostVector());
				LorentzTransform to_H2_CoM(-q2_lab.boostVector());

				FourMomentum q11_H1 = to_H1_CoM.transform(q11_lab);
				FourMomentum q21_H2 = to_H2_CoM.transform(q21_lab);

				cos_theta1 = q1.p3().dot(q11_H1.p3()) / sqrt(q1.p3().mod2()*q11_H1.p3().mod2());
				cos_theta2 = q2.p3().dot(q21_H2.p3()) / sqrt(q2.p3().mod2()*q21_H2.p3().mod2());
			}
			//@}


		private:

			// Data members like post-cuts event weight counters go here


		private:
			TFile* _fOut;
			TTree* _tree;
			double _pt12, _pt34, _eta12, _eta34, _phi12, _phi34;// _dr12, _dr34;
			double _ptX, _etaX, _phiX, _yX;
			int _flav1, _flav2, _flav3, _flav4;
			double _pt1, _eta1, _m1, _phi1;
			double _pt2, _eta2, _m2, _phi2;
			double _pt3, _eta3, _m3, _phi3;
			double _pt4, _eta4, _m4, _phi4;
			double _genWeight, _btagWeight, _weight, _mX, _m12, _m34, _cosThetaStar, _Phi, _Phi1, _cosTheta1, _cosTheta2;
			double _mW12, _mW34, _mt12, _mt34, _dRW12, _dRW34;
			/// @name Histograms
			//@{
			map<string, shared_ptr<YODA::Histo1D> > _histograms_1d;
			map<string, shared_ptr<YODA::Histo2D> > _histograms_2d;
			//@}


	};



	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(Summer_2014_Study);

}
