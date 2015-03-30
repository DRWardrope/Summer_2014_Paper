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
//				std::cout<<"Handler name = "<< handler().runName() << std::endl;	
				//setNeedsCrossSection();
			}
        //int nEventsAccepted = 0;
        
        struct quark_with_dR {
            double quark;
            double index;
            double dR;
        };

		public:

			/// @name Analysis methods
			//@{

			/// Book histograms and initialise projections before the run
			void init()
			{
				//_crossSection = handler().crossSection();
				FinalStateWithGhosts fs(-5.0, 5.0);
				fs.ghostIdPair(PID::BQUARK);
				fs.ghostIdPair(PID::TAU);
				addProjection(fs, "FS");

				FastJets akt4(fs, FastJets::ANTIKT, 0.4);
				//fj.useInvisibles();
				addProjection(akt4, "AntiKt4");

				FastJets ca12(fs, FastJets::CAM, 1.2);
				addProjection(ca12, "CA12");
				bookTree();
				/// @todo Initialise and register projections here

				/// @todo Book histograms here, e.g.:
				// _h_XXXX = bookProfile1D(1, 1, 1);
				// _h_YYYY = bookHisto1D(2, 1, 1);
                
                pidFinder tqf(-5.0, 5.0);
                tqf.ghostIdPair(PID::TQUARK);
                addProjection(tqf, "tqf");

			}


			/// Perform the per-event analysis
			void analyze(const Event& event) 
			{
				const double weight = event.weight();
				// Signal
			/*	ParticleVector higgses, bquarks;
				HepMC::GenEvent::particle_const_iterator it;
				for (it=event.genEvent()->particles_begin(); it!=event.genEvent()->particles_end(); ++it) 
				{

					const HepMC::GenParticle *particle = *it;
					if (particle->pdg_id() == 25) {

						HepMC::GenVertex *decay = particle->end_vertex();
						if (decay && decay->particles_out_size() == 2) {
							higgses.push_back(*particle); // implicit conversion HepMC::GenParticle -> Rivet::Particle

							HepMC::GenVertex::particles_out_const_iterator itv;
							for (itv=decay->particles_out_const_begin(); itv!=decay->particles_out_const_end(); ++itv) {
								if (abs((*itv)->pdg_id()) == 5) {
									bquarks.push_back(*findLastInChain(*itv));
								}
							}
						}
					}
				}

				// For background select based on the final state
				// Find 4 b-quarks that are created from something non-b
				// Pair them into two Higgs bosons with sufficient boost (mass?)
				// Trace them to the final b-quarks

				// QCD background?
				if (higgses.size() == 0) {

					ParticleVector initial_b, initial_bbar;

					HepMC::GenEvent::vertex_const_iterator it_vtx;
					for (it_vtx=event.genEvent()->vertices_begin(); it_vtx!=event.genEvent()->vertices_end(); ++it_vtx) {
						const HepMC::GenVertex *vertex  =*it_vtx;

						bool skip = false;
						HepMC::GenVertex::particles_in_const_iterator it_in;
						for (it_in=vertex->particles_in_const_begin(); it_in!=vertex->particles_in_const_end(); ++it_in) {
							const HepMC::GenParticle *particle = *it_in;

							if (abs(particle->pdg_id()) == 5) skip = true;
						}

						if (skip) continue;

						HepMC::GenVertex::particles_out_const_iterator it_out;
						for (it_out=vertex->particles_out_const_begin(); it_out!=vertex->particles_out_const_end(); ++it_out) {
							const HepMC::GenParticle *particle = *it_out;

							if (particle->pdg_id() == 5) {
								initial_b.push_back(*particle);
							} else if (particle->pdg_id() == -5) {
								initial_bbar.push_back(*particle);
							}
						}
					}

					if (initial_b.size() == 2 && initial_bbar.size() == 2) {
						// Alternative combination - have to check it is the same
						bool alt = false;

						higgses.push_back(Particle(0, initial_b[0].momentum()+initial_bbar[alt ? 1 : 0].momentum()));
						higgses.push_back(Particle(0, initial_b[1].momentum()+initial_bbar[alt ? 0 : 1].momentum()));

						bquarks.push_back(*findLastInChain(initial_b[0].genParticle()));
						bquarks.push_back(*findLastInChain(initial_bbar[alt ? 1 : 0].genParticle()));
						bquarks.push_back(*findLastInChain(initial_b[1].genParticle()));
						bquarks.push_back(*findLastInChain(initial_bbar[alt ? 0 : 1].genParticle()));
					}
				}

				const Jets &akt4 = applyProjection<JetAlg>(event, "AntiKt4").jetsByPt(40*GeV);
				const Jets &ca12 = applyProjection<JetAlg>(event, "CA12"   ).jetsByPt(80*GeV);

				if (higgses.size() == 2 && bquarks.size() == 4 
						&& higgses[0].momentum().pT() > 150*GeV && higgses[1].momentum().pT() > 150*GeV) {

					for (unsigned int i=0; i<4; ++i) {
						fillKinematics("parton", "Parton", bquarks[i].momentum(), weight);
					}
					fillMelaAngles("parton", "Parton", higgses[0].momentum(), bquarks[0].momentum(),  bquarks[1].momentum(), 
							higgses[1].momentum(), bquarks[2].momentum(), bquarks[3].momentum(), weight);

					//PRINT(higgses.size());
					//PRINT(bquarks.size());

					fillReconstructionEfficiency("akt4", "Anti-$k_t$ 0.4 reconstruction efficiency",       higgses[0], bquarks[0], bquarks[1], akt4, weight);
					fillReconstructionEfficiency("akt4", "Anti-$k_t$ 0.4 reconstruction efficiency",       higgses[1], bquarks[2], bquarks[3], akt4, weight);
					fillReconstructionEfficiency("ca12", "Cambridge-Aachen 1.2 reconstruction efficiency", higgses[0], bquarks[0], bquarks[1], ca12, weight);
					fillReconstructionEfficiency("ca12", "Cambridge-Aachen 1.0 reconstruction efficiency", higgses[1], bquarks[2], bquarks[3], ca12, weight);
				}*/
				bookAndFill1D("CutFlow", "Yields vs Selection Step", 0, "Selection Step", 5, -0.5, 4.5, weight);
				const Jets &akt4 = applyProjection<JetAlg>(event, "AntiKt4").jetsByPt(20*GeV); //Lowered from 40 for ttbar veto tests
				std::vector<TagJet> bjets;
				bjets.reserve(4);
				foreach(const Jet &j, akt4) {
					if (j.momentum().pT() > 40*GeV 
							&& fabs(j.momentum().rapidity()) < 2.5)
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
                    /*std::cout << "Event Accepted" << std::endl;
                    std::cout << "Number of jets found = " << bjets.size() << std::endl;
                    std::cout << "Top 4 jets:" << std::endl;
                    for (int j = 0; j < 4; ++j) {
                        std::cout << j <<": pt = "<< bjets[j].pt() <<" GeV, eta = "<< bjets[j].eta() <<", phi = "<< bjets[j].phi() << std::endl;
                        //std::cout << " Momentum: " << bjets[j].p3() << std::endl;
                    }*/
                    
					//Xtt investigation
                    ParticleVector tDaughters;
                    ParticleVector wDaughters;
					//if(fabs(dijet1.mass() - 115.) < 25. && fabs(dijet2.mass() - 110.) < 25.)
					//{
						HepMC::GenEvent::particle_const_iterator it;
						for (it=event.genEvent()->particles_begin(); it!=event.genEvent()->particles_end(); ++it) 
						{
							const HepMC::GenParticle *particle = *it;
							if (abs(particle->pdg_id()) == 24)
							{
							/*	HepMC::GenVertex *decay = particle->end_vertex();
								if (decay && decay->particles_out_size() == 2) 
								{
								//	wBosons.push_back(*particle);
									HepMC::GenVertex::particles_out_const_iterator itv;
									ParticleVector wDaughters;
									for (itv=decay->particles_out_const_begin(); itv!=decay->particles_out_const_end(); ++itv) {
										if (abs((*itv)->pdg_id()) == 5) {
											wDaughters.push_back(*findLastInChain(*itv));
										}
									}
									std::cout<<"Number of w daughters found = "<< wDaughters.size() << std::endl;
									for(int w = 0; w < wDaughters.size(); ++w)
									{
										std::cout<< w <<": id = "<< wDaughters[w].pdgId() <<", pt = "<< wDaughters[w].pt() <<" GeV, eta = "<< wDaughters[w].eta()
												<<", phi = "<< wDaughters[w].phi() << std::endl;
									}
									//std::cout<<"DeltaR(wDaughters[0], wDaughters[1]) = "<< deltaR(wDaughters[0], wDaughters[1]) << std::endl;
									//bookAndFill1D("Wdaughters_dR", "\\Delta R(W daughter 1, W daughter 2)", deltaR(wDaughters[0], wDaughters[1]) ,   "\\Delta R",   50, 0.,   5.,   btaggedWeight);
								}*/
							}else if(abs(particle->pdg_id()) == 6)
							{
								HepMC::GenVertex *decay = particle->end_vertex();
								if (decay && decay->particles_out_size() == 2) 
								{
							//		wBosons.push_back(*particle);
                                    //std::cout << "This top decays!" << std::endl;
                                    Particle p(*particle);
                                    //std::cout << "ID: " << p.pid() << " mass: " << p.mass() << " GeV, pt: " << p.pt() << " GeV, eta: " << p.eta() << ", phi: " << p.phi() << std::endl;
                                    //std::cout << "Momentum: " << p.p3() << std::endl;
									HepMC::GenVertex::particles_out_const_iterator itv;
									//ParticleVector bQuarks;
                                    //ParticleVector wBosons;
                                    //std::cout << "Decay Products:" << std::endl;
									for (itv=decay->particles_out_const_begin(); itv!=decay->particles_out_const_end(); ++itv) {
										//const HepMC::GenParticle *daughter = *itv;
                                        //Particle d(*daughter);
                                        //std::cout << "ID: " << d.pid() << " mass: " << d.mass() << " GeV, pt: " << d.pt() << " GeV, eta: " << d.eta() << ", phi: " << d.phi() << std::endl;
                                        //std::cout << "Momentum: " << d.p3() << std::endl;
                                        if (abs((*itv)->pdg_id()) == 5) {
											//bQuarks.push_back(*itv);
                                            tDaughters.push_back(*itv);
										}
                                        if (abs((*itv)->pdg_id()) == 24) {
                                            const HepMC::GenParticle *wBoson = findLastInChain(*itv);
                                            //wBosons.push_back(*wBoson);
                                            HepMC::GenVertex *decayV = wBoson->end_vertex();
                                            if (decayV && decayV->particles_out_size() == 2)
                                            {
                                                //		wBosons.push_back(*particle);
                                                //std::cout << "This W decays!" << std::endl;
                                                HepMC::GenVertex::particles_out_const_iterator itv2;
                                                for (itv2=decayV->particles_out_const_begin(); itv2!=decayV->particles_out_const_end(); ++itv2) {
                                                    const HepMC::GenParticle *wDaughter = *itv2;
                                                    //const HepMC::GenParticle *wDaughter_debug = findLastInChain(*itv2);
                                                    //Particle wd(*wDaughter);
                                                    wDaughters.push_back(*wDaughter);
                                                    if (abs((*itv2)->pdg_id()) <= 6) { //if it is a hadron, otherwise do not bother matching
                                                        tDaughters.push_back(*wDaughter);
                                                    }
                                                    //std::cout << "  ID: " << wd.pid() << " mass: " << wd.mass() << " GeV, pt: " << wd.pt() << " GeV, eta: " << wd.eta() << ", phi: " << wd.phi() << std::endl;
                                                    //std::cout << "  Momentum: " << wd.p3() << std::endl;
                                                }
                                            }
                                        }
									}
                                    
                                    /*
									std::cout<<"Number of b quarks found = "<< bQuarks.size() << std::endl;
									for(int w = 0; w < bQuarks.size(); ++w)
									{
										std::cout<< w <<": id = "<< bQuarks[w].pdgId() <<", pt = "<< bQuarks[w].pt() <<" GeV, eta = "<< bQuarks[w].eta()
												<<", phi = "<< bQuarks[w].phi() << std::endl;
									}
                                    
                                    std::cout<<"Number of w Bosons found = "<< wBosons.size() << std::endl;
                                    for(int w = 0; w < wBosons.size(); ++w)
                                    {
                                        std::cout<< w <<": id = "<< wBosons[w].pdgId() <<", pt = "<< wBosons[w].pt() <<" GeV, eta = "<< wBosons[w].eta()
                                        <<", phi = "<< wBosons[w].phi() << std::endl;
                                    }
                                    
                                    std::cout<<"Number of w daughters found = "<< wDaughters.size() << std::endl;
                                    for(int w = 0; w < wDaughters.size(); ++w)
                                    {
                                        std::cout<< w <<": id = "<< wDaughters[w].pdgId() <<", pt = "<< wDaughters[w].pt() <<" GeV, eta = "<< wDaughters[w].eta()
                                        <<", phi = "<< wDaughters[w].phi() << std::endl;
                                    }
                                    */
									//std::cout<<"DeltaR(tDaughters[0], tDaughters[1]) = "<< deltaR(tDaughters[0], tDaughters[1]) << std::endl;
									//bookAndFill1D("tDaughters_dR", "\\Delta R(t daughter 1, t daughter 2)", deltaR(tDaughters[0], tDaughters[1]) ,   "\\Delta R",   50, 0.,   5.,   btaggedWeight);
								}

							}
						}
					//}
                    //Now that the event has passed the dijet selection, extract all top quarks from event and add properties of last pair to tree
                    const FinalState& tqf = applyProjection<FinalState>(event, "tqf");
                    fillPidVar(tqf, PID::TQUARK);
                    /*
                    std::cout<<"Number of t daughter hadrons found = "<< tDaughters.size() << std::endl;
                    for(int w = 0; w < tDaughters.size(); ++w)
                    {
                        std::cout<< w <<": id = "<< tDaughters[w].pdgId() <<", pt = "<< tDaughters[w].pt() <<" GeV, eta = "<< tDaughters[w].eta()
                        <<", phi = "<< tDaughters[w].phi() << std::endl;
                    }
                    */
                    bJetMatchPar(bjets, tDaughters, wDaughters);
                    
					fillKinematics("dijet", "Dijet", dijet1, btaggedWeight);
					fillKinematics("dijet", "Dijet", dijet2, btaggedWeight);
					for (unsigned int i=0; i<4; ++i) {
						fillKinematics("jet", "Jet", bjets[i].momentum(), btaggedWeight);
					}
					fillXttVariables(akt4,"jet", "Jet", dijet1, bjets[0], bjets[1], dijet2, bjets[2], bjets[3], btaggedWeight, tDaughters);
					fillMelaAngles(event, "jet", "Jet", dijet1, bjets[0], bjets[1], dijet2, bjets[2], bjets[3], btaggedWeight);
					//dijetSelection only satisfied if ≥2 jets are found. Use weight, since x=2 represents 2 dijets w/o b-tagging req.
					bookAndFill1D("CutFlow", "Yields vs Selection Step", 2, "Selection Step", 5, -0.5, 4.5, weight);
					//Use btaggedWeight to replicate effect of b-tagging.
					bookAndFill1D("CutFlow", "Yields vs Selection Step", 3, "Selection Step", 5, -0.5, 4.5, btaggedWeight);
                    
                    //nEventsAccepted++;
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
                    const YODA::Histo1D Hist0 = *hIt->second;
					const YODA::HistoBin1D& bin0 = Hist0.binAt(0);
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
                
                //std::cout << "Number of events Accepted: " << nEventsAccepted << std::endl;

			}
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
                //For the TQUARKS: A is t, B is tbar
                _tree->Branch("ptA", &_ptA, "ptA/D");
                _tree->Branch("mA", &_mA, "mA/D");
                _tree->Branch("etaA", &_etaA, "etaA/D");
                _tree->Branch("yA", &_yA, "yA/D");
                _tree->Branch("phiA", &_phiA, "phiA/D");
                
                _tree->Branch("ptB", &_ptB, "ptB/D");
                _tree->Branch("mB", &_mB, "mB/D");
                _tree->Branch("etaB", &_etaB, "etaB/D");
                _tree->Branch("yB", &_yB, "yB/D");
                _tree->Branch("phiB", &_phiB, "phiB/D");
                //Jet Matching related
                _tree->Branch("quarkMiss12", &_quarkMiss12, "quarkMiss12/D");
                _tree->Branch("ptMiss12", &_ptMiss12, "ptMiss12/D");
                _tree->Branch("mMiss12", &_mMiss12, "mMiss12/D");
                _tree->Branch("etaMiss12", &_etaMiss12, "etaMiss12/D");
                _tree->Branch("yMiss12", &_yMiss12, "yMiss12/D");
                _tree->Branch("phiMiss12", &_phiMiss12, "phiMiss12/D");
                
                _tree->Branch("quarkMiss34", &_quarkMiss34, "quarkMiss34/D");
                _tree->Branch("ptMiss34", &_ptMiss34, "ptMiss34/D");
                _tree->Branch("mMiss34", &_mMiss34, "mMiss34/D");
                _tree->Branch("etaMiss34", &_etaMiss34, "etaMiss34/D");
                _tree->Branch("yMiss34", &_yMiss34, "yMiss34/D");
                _tree->Branch("phiMiss34", &_phiMiss34, "phiMiss34/D");
                
                _tree->Branch("charge12", &_charge12, "charge12/D");
                _tree->Branch("charge34", &_charge34, "charge34/D");
                
                _tree->Branch("bLeadStatus12", &_bLeadStatus12, "bLeadStatus12/D");
                _tree->Branch("bLeadStatus34", &_bLeadStatus34, "bLeadStatus34/D");
                _tree->Branch("bothBleading", &_bothBleading, "bothBleading/D");
                _tree->Branch("dijetsHaveTwoBquarks", &_dijetsHaveTwoBquarks, "dijetsHaveTwoBquarks/D");
                
                _tree->Branch("dR_sublead12", &_dR_sublead12, "dR_sublead12/D");
                _tree->Branch("dR_sublead34", &_dR_sublead34, "dR_sublead34/D");
                _tree->Branch("dR_lead12", &_dR_lead12, "dR_lead12/D");
                _tree->Branch("dR_lead34", &_dR_lead34, "dR_lead34/D");
                
                _tree->Branch("dR_Wquarks12", &_dR_Wquarks12, "dR_Wquarks12/D");
                _tree->Branch("dR_Wquarks34", &_dR_Wquarks34, "dR_Wquarks34/D");
                
                _tree->Branch("numJetsMatched", &_numJetsMatched, "numJetsMatched/D");
                _tree->Branch("numQuarks", &_numQuarks, "numQuarks/D");
                _tree->Branch("quarkMatchStatus", &_quarkMatchStatus, "quarkMatchStatus/D");
                _tree->Branch("alignStatus", &_alignStatus, "alignStatus/D");
                _tree->Branch("multiMatch", &_multiMatch, "multiMatch/D");
                
                _tree->Branch("quark1", &_quark1, "quark1/D");
                _tree->Branch("quark2", &_quark2, "quark2/D");
                _tree->Branch("quark3", &_quark3, "quark3/D");
                _tree->Branch("quark4", &_quark4, "quark4/D");
                
                _tree->Branch("deltaR1", &_deltaR1, "deltaR1/D");
                _tree->Branch("deltaR2", &_deltaR2, "deltaR2/D");
                _tree->Branch("deltaR3", &_deltaR3, "deltaR3/D");
                _tree->Branch("deltaR4", &_deltaR4, "deltaR4/D");
                //3rd Jet related
                _tree->Branch("numNonDup", &_numNonDup, "numNonDup/D");
                _tree->Branch("ptMinAll", &_ptMinAll, "ptMinAll/D");
                _tree->Branch("etaMaxAll", &_etaMaxAll, "etaMaxAll/D");
                _tree->Branch("numClosestFound", &_numClosestFound, "numClosestFound/D");
                _tree->Branch("multiMatch3rd", &_multiMatch3rd, "multiMatch3rd/D");
                _tree->Branch("num3rdJetsMatched", &_num3rdJetsMatched, "num3rdJetsMatched/D");
                _tree->Branch("sameClosest", &_sameClosest, "sameClosest/D");
                _tree->Branch("3rdMatchStatus", &_3rdMatchStatus, "3rdMatchStatus/D");
                
                _tree->Branch("quarkMatched3rdJet12", &_quarkMatched3rdJet12, "quarkMatched3rdJet12/D");
                _tree->Branch("pt3rdJet12", &_pt3rdJet12, "pt3rdJet12/D");
                _tree->Branch("m3rdJet12", &_m3rdJet12, "m3rdJet12/D");
                _tree->Branch("eta3rdJet12", &_eta3rdJet12, "eta3rdJet12/D");
                _tree->Branch("y3rdJet12", &_y3rdJet12, "y3rdJet12/D");
                _tree->Branch("phi3rdJet12", &_phi3rdJet12, "phi3rdJet12/D");
                
                _tree->Branch("quarkMatched3rdJet34", &_quarkMatched3rdJet34, "quarkMatched3rdJet34/D");
                _tree->Branch("pt3rdJet34", &_pt3rdJet34, "pt3rdJet34/D");
                _tree->Branch("m3rdJet34", &_m3rdJet34, "m3rdJet34/D");
                _tree->Branch("eta3rdJet34", &_eta3rdJet34, "eta3rdJet34/D");
                _tree->Branch("y3rdJet34", &_y3rdJet34, "y3rdJet34/D");
                _tree->Branch("phi3rdJet34", &_phi3rdJet34, "phi3rdJet34/D");
                
                _tree->Branch("dR3rdJet12", &_dR3rdJet12, "dR3rdJet12/D");
                _tree->Branch("dR3rdJet34", &_dR3rdJet34, "dR3rdJet34/D");
                _tree->Branch("dRMatched3rdJet12", &_dRMatched3rdJet12, "dRMatched3rdJet12/D");
                _tree->Branch("dRMatched3rdJet34", &_dRMatched3rdJet34, "dRMatched3rdJet34/D");
                
                _tree->Branch("numClosestPass12", &_numClosestPass12, "numClosestPass12/D");
                _tree->Branch("numClosestPass34", &_numClosestPass34, "numClosestPass34/D");

			}
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


				/*if (fabs(H1.pT() - 200*GeV) < 10*GeV) {
					const FourMomentum &softest1 = b11.pT() < b12.pT() ? b11 : b12;
					bookAndFill2D(label+"_ptmin_vs_cos_theta", title+" $p_{T,min}$ vs. $|\\cos(\\theta)|$", softest1.pT()/GeV, "$p_{T,min}$", 50, 0., 250., fabs(cos_theta1), "$|\\cos(\\theta)|$", 50, 0., 1., weight);
					bookAndFill2D(label+"_ptmin_vs_dr",        title+" $p_{T,min}$ vs. $\\Delta R(b,b)$",   softest1.pT()/GeV, "$p_{T,min}$", 50, 0., 250., deltaR(b11, b12), "$\\Delta R(b,b)$",   50, 0., 5., weight);
				}

				if (fabs(H2.pT() - 200*GeV) < 10*GeV) {
					const FourMomentum &softest2 = b21.pT() < b22.pT() ? b21 : b22;
					bookAndFill2D(label+"_ptmin_vs_cos_theta", title+" $p_{T,min}$ vs. $|\\cos(\\theta)|$", softest2.pT()/GeV, "$p_{T,min}$", 50, 0., 250., fabs(cos_theta2), "$|\\cos(\\theta)|$", 50, 0., 1., weight);
					bookAndFill2D(label+"_ptmin_vs_dr",        title+" $p_{T,min}$ vs. $\\Delta R(b,b)$",   softest2.pT()/GeV, "$p_{T,min}$", 50, 0., 250., deltaR(b21, b22), "$\\Delta R(b,b)$",   50, 0., 5., weight);
				}*/
			}
            //Based on fillMelaAngles, Fill properties of particle pair given a pid
            void fillPidVar(const FinalState& fs, PdgId pid)
            {
                //A (incomplete) subscript mess that works, keep at the moment in case iterator version goes crazy
                /*
                double i = fs.particles().size()-1;
                bool isTQUARK = false;
                double index = fs.particles()[i].pid();
                if ( index == PID::TQUARK ) {
                    std::cout << "I FOUND A TQUARK FIRST" << std::endl;
                    isTQUARK = true;
                }
                else if ( index == -PID::TQUARK ) {
                    std::cout << "I FOUND A ANTITQUARK FIRST" << std::endl;
                    isTQUARK = true;
                }
                else {
                    std::cout << "ITS NOT A TQUARK" << std::endl;
                }
                if (isTQUARK == true) {
                    std::cout << "Particle number: " << i+1 << std::endl;
                    std::cout << "ID: " << fs.particles()[i].pid() << " mass: " << fs.particles()[i].mass() << " pt: " << fs.particles()[i].pt() << " eta: " << fs.particles()[i].eta() << " phi: " << fs.particles()[i].phi() << std::endl;
                    //std::cout << "testing first" << std::endl;
                }
                 */
                
                //Let's use an iterator to start at the end! ...I wanted to use a reverse iterator but i just couldn't get it to work so here's an interim solution
                //Consider throwing exceptions if there are no particles in projection or particle pair not found? (Since a single instance of either of these can invalidate the data)
                bool filledPar = false;
                bool filledAntiPar = false;
                if (fs.particles().size() == 0) {
                    std::cout << "There are no particles in this projection!" << std::endl;
                    throw 10;
                }
                for (Particles::const_iterator it = fs.particles().end()-1; it!=fs.particles().begin(); --it) {
                    Particle p(*it);
                    if ( p.pid() == pid && filledPar == false ) {
                        //std::cout << "I FOUND A" << std::endl;
                        _ptA = p.pt();
                        _mA = p.mass();
                        _etaA = p.eta();
                        _yA = p.rap();
                        _phiA = p.phi();
                        filledPar = true;
                    }
                    else if ( p.pid() == -pid && filledAntiPar == false) {
                        //std::cout << "I FOUND B(ANTI)" << std::endl;
                        _ptB = p.pt();
                        _mB = p.mass();
                        _etaB = p.eta();
                        _yB = p.rap();
                        _phiB = p.phi();
                        filledAntiPar = true;
                    }
                    if (filledPar && filledAntiPar) {
                        //std::cout << "Particle pair found!" << std::endl;
                        break;
                    }

                }
                if (filledPar == false || filledAntiPar == false) {
                    std::cout << "Particle pair not found..." << std::endl;
                    throw 20;
                }
            }
        
            //Matches particles found from tracing decays to identified jets. Currently optimised for ttbar. Sorry about the magic numbers.
            //Also does diagnostics for 6 hardon 4 jets matched events. Wanted to move it to a seperate function, but using shared variables for now.
            void bJetMatchPar(std::vector<TagJet> &bjets, ParticleVector &tDaughters, ParticleVector &wDaughters)
            {
                
                std::map<int, quark_with_dR> quark;
                double multiMatch = 0;
                double numJetsMatched = 0;
                int numToMatch = 4;
                numJetsMatched = matchJetToPar(numToMatch, multiMatch, quark, bjets, tDaughters);
                _multiMatch = multiMatch;
 
                //Determining Match status code based on number of hardrons found vs number of jets found: Ideally get 3 for majority of events, 1 for rare tau decay
                //2 signifies imperfect matching despite fully hadronic decay: preliminary investigation suggests this is not due to taking first particle in iterator.
                if (numJetsMatched == -5) {
                    _quarkMatchStatus = -5;
                }
                else if (tDaughters.size() == 6 && numJetsMatched == 4) {
                    _quarkMatchStatus = 3;
                }
                else if (tDaughters.size() == 6 && numJetsMatched < 4) {
                    if (tDaughters.size() == 6 && numJetsMatched == 0) {
                        _quarkMatchStatus = -1;
                    }
                    else {
                        _quarkMatchStatus = 2;
                    }
                }
                else if (tDaughters.size() == 4 && numJetsMatched == 4) {
                    _quarkMatchStatus = -2;
                }
                else if (tDaughters.size() == 4 && numJetsMatched < 4) {
                    if (tDaughters.size() == 4 && numJetsMatched == 0) {
                        _quarkMatchStatus = -3;
                    }
                    else {
                       _quarkMatchStatus = 1;
                    }
                }
                else if (tDaughters.size() == 2 && numJetsMatched <= 2) {
                    if (tDaughters.size() == 2 && numJetsMatched == 0) {
                        _quarkMatchStatus = -4;
                    }
                    else {
                        _quarkMatchStatus = 0;
                    }
                }
                else {
                    _quarkMatchStatus = -6;
                }
                //thirdIndex and MatchStatus will be passed to fillXttVariables, hence resetting here
                thirdIndex_A = -1;
                thirdIndex_B = -1;
                //If fully hadronic and perfect match, perform diagnosis
                if (tDaughters.size() == 6 && numJetsMatched == 4) {
                    //Determine dijet charge from associated particles
                    const double upcharge = 2./3;
                    const double downcharge = -1./3;
                    double charge12 = 0;
                    for (int i = 1; i < 3; i++) {
                        double is_negative = 1;
                        if (quark[i].quark < 0) {
                            is_negative = -1;
                        }
                        if (fabs(quark[i].quark) == 5) {
                            charge12 = charge12 + ((is_negative)*(downcharge));
                        }
                        else if (fabs(quark[i].quark) == 4) {
                            charge12 = charge12 + ((is_negative)*(upcharge));
                        }
                        else if (fabs(quark[i].quark) == 3) {
                            charge12 = charge12 + ((is_negative)*(downcharge));
                        }
                        else if (fabs(quark[i].quark) == 2) {
                            charge12 = charge12 + ((is_negative)*(upcharge));
                        }
                        else if (fabs(quark[i].quark) == 1) {
                            charge12 = charge12 + ((is_negative)*(downcharge));
                        }
                    }
                    double charge34 = 0;
                    for (int i = 3; i < 5; i++) {
                        double is_negative = 1;
                        if (quark[i].quark < 0) {
                            is_negative = -1;
                        }
                        if (fabs(quark[i].quark) == 5) {
                            charge34 = charge34 + ((is_negative)*(downcharge));
                        }
                        else if (fabs(quark[i].quark) == 4) {
                            charge34 = charge34 + ((is_negative)*(upcharge));
                        }
                        else if (fabs(quark[i].quark) == 3) {
                            charge34 = charge34 + ((is_negative)*(downcharge));
                        }
                        else if (fabs(quark[i].quark) == 2) {
                            charge34 = charge34 + ((is_negative)*(upcharge));
                        }
                        else if (fabs(quark[i].quark) == 1) {
                            charge34 = charge34 + ((is_negative)*(downcharge));
                        }
                    }
                    _charge12 = charge12;
                    _charge34 = charge34;
                    //Orientation tests: is the dijet composed of particles of only 1 top or both tops?
                    bool from_top_12 = (quark[1].index < 3 && quark[2].index < 3);
                    bool from_antitop_12 = (quark[1].index > 2 && quark[2].index > 2);
                    bool from_top_34 = (quark[3].index < 3 && quark[4].index < 3);
                    bool from_antitop_34 = (quark[3].index > 2 && quark[4].index > 2);
                    if ((from_top_12 && from_antitop_34) || (from_top_34 && from_antitop_12)) { //If each dijet only has 1 contributing top, find missing particle
                        _alignStatus = 2;
                        if (from_top_12) {
                            for (double i=0; i < 3; i++) {
                                if (quark[1].index == i) continue;
                                if (quark[2].index == i) continue;
                                thirdIndex_A = i;
                            }
                            _dR_Wquarks12 = deltaR(wDaughters[0], wDaughters[1]);
                        }
                        else if (from_antitop_12) {
                            for (double i=3; i < 6; i++) {
                                if (quark[1].index == i) continue;
                                if (quark[2].index == i) continue;
                                thirdIndex_A = i;
                            }
                            _dR_Wquarks12 = deltaR(wDaughters[2], wDaughters[3]);
                        }
                        //std::cout << "thirdIndex_A = " << thirdIndex_A << std::endl;
                        if (from_top_34) {
                            for (double i=0; i < 3; i++) {
                                if (quark[3].index == i) continue;
                                if (quark[4].index == i) continue;
                                thirdIndex_B = i;
                            }
                            _dR_Wquarks34 = deltaR(wDaughters[0], wDaughters[1]);
                        }
                        else if (from_antitop_34) {
                            for (double i=3; i < 6; i++) {
                                if (quark[3].index == i) continue;
                                if (quark[4].index == i) continue;
                                thirdIndex_B = i;
                            }
                            _dR_Wquarks34 = deltaR(wDaughters[2], wDaughters[3]);
                        }
                        //std::cout << "thirdIndex_B = " << thirdIndex_B << std::endl;
                        if (thirdIndex_A == -1 || thirdIndex_B == -1) { //If either third particle not identifed for any reason, crash
                            std::cout << "Third particle not found!" << std::endl;
                            throw 40;
                        }
                        //Third particle and lead/subleading jet dR calculations
                        if (bjets[0].tagEff() >= bjets[1].tagEff()) {
                            _dR_sublead12 = deltaR(bjets[1], tDaughters[thirdIndex_A]);
                            _dR_lead12 = deltaR(bjets[0], tDaughters[thirdIndex_A]);
                            if (fabs(quark[1].quark) == 5) {
                                _bLeadStatus12 = 2;
                            }
                            else if (fabs(quark[2].quark) == 5) {
                                _bLeadStatus12 = 1;
                            }
                            else {
                                _bLeadStatus12 = 0;
                            }
                        }
                        else {
                            _dR_sublead12 = deltaR(bjets[0], tDaughters[thirdIndex_A]);
                            _dR_lead12 = deltaR(bjets[1], tDaughters[thirdIndex_A]);
                            if (fabs(quark[2].quark) == 5) {
                                _bLeadStatus12 = 2;
                            }
                            else if (fabs(quark[1].quark) == 5) {
                                _bLeadStatus12 = 1;
                            }
                            else {
                                _bLeadStatus12 = 0;
                            }
                        }
                        
                        if (bjets[2].tagEff() >= bjets[3].tagEff()) {
                            _dR_sublead34 = deltaR(bjets[3], tDaughters[thirdIndex_B]);
                            _dR_lead34 = deltaR(bjets[2], tDaughters[thirdIndex_B]);
                            if (fabs(quark[3].quark) == 5) {
                                _bLeadStatus34 = 2;
                            }
                            else if (fabs(quark[4].quark) == 5) {
                                _bLeadStatus34 = 1;
                            }
                            else {
                                _bLeadStatus34 = 0;
                            }
                        }
                        else {
                            _dR_sublead34 = deltaR(bjets[2], tDaughters[thirdIndex_B]);
                            _dR_lead34 = deltaR(bjets[3], tDaughters[thirdIndex_B]);
                            if (fabs(quark[4].quark) == 5) {
                                _bLeadStatus34 = 2;
                            }
                            else if (fabs(quark[3].quark) == 5) {
                                _bLeadStatus34 = 1;
                            }
                            else {
                                _bLeadStatus34 = 0;
                            }
                        }
                        /*
                        std::cout << "bLeadStatus12: " << _bLeadStatus12 << std::endl;
                        std::cout << "bLeadStatus34: " << _bLeadStatus34 << std::endl;
                        std::cout << "dRWquarks12: " << _dR_Wquarks12 << std::endl;
                        std::cout << "dRWquarks34: " << _dR_Wquarks34 << std::endl;
                         */
                        if (_bLeadStatus12 == 2 && _bLeadStatus34 == 2) {
                            _bothBleading = 2;
                        }
                        else if (_bLeadStatus12 == 2 || _bLeadStatus34 == 2) {
                            _bothBleading = 1;
                        }
                        else {
                            _bothBleading = 0;
                        }
                        if (_bLeadStatus12 > 0 && _bLeadStatus34 > 0) {
                            _dijetsHaveTwoBquarks = 2;
                        }
                        else if (_bLeadStatus12 > 0 || _bLeadStatus34 > 0) {
                            _dijetsHaveTwoBquarks = 1;
                        }
                        else {
                            _dijetsHaveTwoBquarks = 0;
                        }
                        //Third particle information storage:
                        _quarkMiss12 = tDaughters[thirdIndex_A].pdgId();
                        _ptMiss12 = tDaughters[thirdIndex_A].pt();
                        _etaMiss12 = tDaughters[thirdIndex_A].eta();
                        _phiMiss12 = tDaughters[thirdIndex_A].phi();
                        _yMiss12 = tDaughters[thirdIndex_A].rap();
                        _mMiss12 = tDaughters[thirdIndex_A].mass();
                        
                        _quarkMiss34 = tDaughters[thirdIndex_B].pdgId();
                        _ptMiss34 = tDaughters[thirdIndex_B].pt();
                        _etaMiss34 = tDaughters[thirdIndex_B].eta();
                        _phiMiss34 = tDaughters[thirdIndex_B].phi();
                        _yMiss34 = tDaughters[thirdIndex_B].rap();
                        _mMiss34 = tDaughters[thirdIndex_B].mass();
                    }
                    else if (from_top_12 || from_antitop_34 || from_top_34 || from_antitop_12) { //If only 1 dijet is compeletely from 1 top
                        _alignStatus = 1;
                    }
                    else {
                        _alignStatus = 0;
                    }
                }
                else { //Fill unused variables with irrelevant constants
                    _charge12 = -5;
                    _charge34 = -5;
                    _alignStatus = -1;
                    _dR_sublead12 = -10;
                    _dR_sublead34 = -10;
                    _dR_lead12 = -10;
                    _dR_lead34 = -10;
                    _bLeadStatus12 = -5;
                    _bLeadStatus34 = -5;
                    _bothBleading = -5;
                    _dijetsHaveTwoBquarks = -5;
                    _dR_Wquarks12 = -10;
                    _dR_Wquarks34 = -10;
                    
                    _quarkMiss12 = -10;
                    _ptMiss12 = -99;
                    _etaMiss12 = -99;
                    _phiMiss12 = -99;
                    _yMiss12 = -99;
                    _mMiss12 = -99;
                    
                    _quarkMiss34 = -10;
                    _ptMiss34 = -99;
                    _etaMiss34 = -99;
                    _phiMiss34 = -99;
                    _yMiss34 = -99;
                    _mMiss34 = -99;
                }
                
                _numQuarks = tDaughters.size();
                _numJetsMatched = numJetsMatched;
                
                //quark1 = fabs(quark[1].quark);
                _quark1 = quark[1].quark;
                _deltaR1 = quark[1].dR;
                //_quark2 = fabs(quark[2].quark);
                _quark2 = quark[2].quark;
                _deltaR2 = quark[2].dR;
                //_quark3 = fabs(quark[3].quark);
                _quark3 = quark[3].quark;
                _deltaR3 = quark[3].dR;
                //_quark4 = fabs(quark[4].quark);
                _quark4 = quark[4].quark;
                _deltaR4 = quark[4].dR;
            }

            double matchJetToPar(int numToMatch, double &multiMatch, std::map<int, quark_with_dR> &quark, std::vector<TagJet> &bjets, ParticleVector &tDaughters)
            {
                double numJetsMatched = 0;
                std::vector<int> barcodes;
                int barCount = 0;
                quark_with_dR nullQuark;
                nullQuark.quark = -10;
                nullQuark.index = -1;
                nullQuark.dR = -10;
                try {
                    for (int i=0; i<numToMatch; ++i) {
                        std::map<int, double> dR;
                        quark_with_dR matchingQuark;
                        int minkey = -1;
                        for (int j=0; j<tDaughters.size(); ++j) {
                            //std::cout << "dR for jet no. " << i << " and tDaughter no. " << j << " = " << deltaR(bjets[i], tDaughters[j]) << std::endl;
                            double dR_temp = deltaR(bjets[i], tDaughters[j]);
                            if (dR_temp < 0.3) {
                                //std::cout << "dR for jet no. " << i << " and tDaughter no. " << j << " = " << deltaR(bjets[i], tDaughters[j]) << std::endl;
                                dR[j] = (dR_temp);
                                if (minkey == -1) {
                                    minkey = j;
                                } else if (dR_temp < dR[minkey]){
                                    //std::cout << "Found a smaller dR!" << std::endl;
                                    minkey = j;
                                }
                            }
                            
                        }
                        /*if (dR.size() == 0) {
                         std::cout << "No match for jet no. " << i << std::endl;
                         }
                         if (dR.size() == 1) {
                         //std::cout << "Perfect!" << std::endl;
                         }
                         if (dR.size() >= 1) {
                         std::cout << "Matched!" << std::endl;
                         }*/
                        if (dR.size() > 1) {
                            //std::cout << "oh noes multiple matches!" << std::endl;
                            multiMatch = 1;
                        }
                        if (minkey == -1) {
                            matchingQuark = nullQuark;
                            quark[i+1] = matchingQuark;
                            //std::cout << "no match" << std::endl;
                            //std::cout << "dR for jet no. " << i << " and tDaughter no. " << minkey << " = " << matchingQuark.dR << ", flavour = " << matchingQuark.quark <<std::endl;
                            continue;
                        }
                        matchingQuark.quark = tDaughters[minkey].pdgId();
                        for (std::vector<int>::iterator it = barcodes.begin() ; it != barcodes.end(); ++it) { //If same particle used to match more than 1 jet, crash!
                            //std::cout << barCount << " " << *it << std::endl;
                            if (tDaughters[minkey].genParticle()->barcode() == *it) {
                                //std::cout << "This particle has been used!" << std::endl;
                                //matchingQuark.dR = dR[minkey];
                                //matchingQuark.index = minkey;
                                //std::cout << "dR for jet no. " << i << " and tDaughter no. " << minkey << " = " << matchingQuark.dR << ", flavour = " << matchingQuark.quark <<std::endl;
                                throw 30;
                            }
                            barCount++;
                        }
                        barcodes.push_back(tDaughters[minkey].genParticle()->barcode());
                        matchingQuark.dR = dR[minkey];
                        matchingQuark.index = minkey;
                        //std::cout << "match" << std::endl;
                        //std::cout << "dR for jet no. " << i << " and tDaughter no. " << minkey << " = " << matchingQuark.dR << ", flavour = " << matchingQuark.quark <<std::endl;
                        //std::cout << tDaughters[minkey].genParticle()->barcode() << std::endl;
                        //std::cout << "--------" << std::endl;
                        quark[i+1] = matchingQuark;
                        numJetsMatched++;
                    }
                }
                catch (int e) {
                    if (e == 30) {
                        //std::cout << "Invalidating event: set to -5 jets matched" << std::endl;
                        numJetsMatched = -5;
                        quark.clear();
                        for (int i=0; i<numToMatch; ++i) {
                            quark[i+1] = nullQuark;
                        }
                    }
                    else {
                        std::cout << "Something went wrong" << std::endl;
                        throw e;
                    }
                }
                /*std::cout << "Results" << std::endl;
                std::cout << "num Matched: " << numJetsMatched << std::endl;
                for (int i=1; i<=numToMatch; ++i) {
                    std::cout << "dR for jet no. " << i-1 << " and tDaughter no. " << quark[i].index << " = " << quark[i].dR << ", flavour = " << quark[i].quark <<std::endl;
                }*/
                return numJetsMatched;
            }
        
            /// Modified fillXttVariables 1/2/2015 to diagonise 3rd jet
            void fillXttVariables(const Jets& jets, std::string label, std::string title, const FourMomentum &H1, const TagJet &b11,
                                  const TagJet &b12, const FourMomentum &H2, const TagJet &b21, const TagJet &b22, double weight, ParticleVector &tDaughters)
            {
                _mW12 = -99, _mW34 = -99, _mt12 = -99, _mt34 = -99, _dRW12 = -99, _dRW34 = -99;
                
                _pt3rdJet12 = -99;
                _eta3rdJet12 = -99;
                _phi3rdJet12 = -99;
                _y3rdJet12 = -99;
                _m3rdJet12 = -99;
                _pt3rdJet34 = -99;
                _eta3rdJet34 = -99;
                _phi3rdJet34 = -99;
                _y3rdJet34 = -99;
                _m3rdJet34 = -99;
                
                _numClosestFound = -10;
                _multiMatch3rd = -10;
                _num3rdJetsMatched = -10;
                _sameClosest = -1;
                _3rdMatchStatus = -4;
                
                _dR3rdJet12 = -99;
                _dR3rdJet34 = -99;
                
                _quarkMatched3rdJet12 = -10;
                _dRMatched3rdJet12 = -99;
                _quarkMatched3rdJet34 = -10;
                _dRMatched3rdJet34 = -99;

                
                Jets noDups;
                int counter = 0;
                double ptMinAll = -99;
                double etaMaxAll = -99;
                foreach(const Jet& j, jets)
                {
                    if (counter == 0) {
                        ptMinAll = j.pt();
                        etaMaxAll = fabs(j.eta());
                    }
                    if (j.pt() < ptMinAll) {
                        ptMinAll = j.pt();
                    }
                    if (fabs(j.eta()) > etaMaxAll) {
                        etaMaxAll = fabs(j.eta());
                    }
                    //std::cout << "A jet: " << counter << std::endl;
                    //std::cout <<" pt = "<< j.pt() <<" GeV, eta = "<< j.eta() <<", phi = "<< j.phi() << std::endl;
                    counter++;
                    if(j.eta() == b11.eta() && j.phi() == b11.phi()) continue;
                    if(j.eta() == b12.eta() && j.phi() == b12.phi()) continue;
                    if(j.eta() == b21.eta() && j.phi() == b21.phi()) continue;
                    if(j.eta() == b22.eta() && j.phi() == b22.phi()) continue;
                    noDups.push_back(j);
                }
                //std::cout << "Min pt: " << ptMinAll << std::endl;
                //std::cout << "Max eta: " << etaMaxAll << std::endl;
                _ptMinAll = ptMinAll;
                _etaMaxAll = etaMaxAll;
                //std::cout << "In Xtt, Match Status: " << _quarkMatchStatus << std::endl;
                //std::cout << "thirdIndexA: " << thirdIndex_A << std::endl;
                //std::cout << "thirdIndexA: " << thirdIndex_B << std::endl;
                assert(noDups.size() == jets.size() - 4);
                const TagJet &leastTagged12 = b11.tagEff() >= b12.tagEff() ? b12 : b11;
                const TagJet &leastTagged34 = b21.tagEff() >= b22.tagEff() ? b22 : b21;
                //float mindR12 = 1.5, mindR34 = 1.5;
                float mindR12 = 2.0, mindR34 = 2.0;
                float mindR12_debug = 2.0, mindR34_debug = 2.0;
                //float mindR12_to1 = 100, mindR34_to1 = 100;
                //float dR12_to1 = 200, dR34_to1 = 200;
                const Jet *closestJet12=NULL, *closestJet34=NULL;
                //std::cout << "Number of non-selected jets found = " << noDups.size() << std::endl;
                _numNonDup = noDups.size();
                std::vector<TagJet> closestJets;
                double numClosestPass12 = 0;
                double numClosestPass34 = 0;
                /*
                foreach(const Jet& j, noDups) //Now using additional condition that dR is closest to 1
                {
                    float dR12 = deltaR(leastTagged12, j);
                    if(dR12 < mindR12 && dR12 > 0.1)
                    {
                        //mindR12 = dR12;
                        dR12_to1 = fabs(1-dR12);
                        if (dR12_to1 < mindR12_to1) {
                            mindR12_to1 = dR12_to1;
                            closestJet12 = &j;
                            mindR12_debug = dR12;
                        }
                        //closestJet12 = &j;
                        numClosestPass12++;
                    }
                    float dR34 = deltaR(leastTagged34, j);
                    if(dR34 < mindR34 && dR34 > 0.1)
                    {
                        //mindR34 = dR34;
                        dR34_to1 = fabs(1-dR34);
                        if (dR34_to1 < mindR34_to1) {
                            mindR34_to1 = dR34_to1;
                            closestJet34 = &j;
                            mindR34_debug = dR34;
                        }
                        //closestJet34 = &j;
                        numClosestPass34++;
                    }
                    //std::cout << "Non dup jet:" << std::endl;
                    //std::cout <<" pt = "<< j.pt() <<" GeV, eta = "<< j.eta() <<", phi = "<< j.phi() << std::endl;
                }
                 */
                //Old use minimun dR for 3 jet code
                foreach(const Jet& j, noDups)
                {
                    float dR12 = deltaR(leastTagged12, j);
                    if(dR12 < mindR12_debug && dR12 > 0.1)
                    {
                        numClosestPass12++;
                    }
                    if(dR12 < mindR12 && dR12 > 0.1)
                    {
                        mindR12 = dR12;
                        closestJet12 = &j;
                    }
                    float dR34 = deltaR(leastTagged34, j);
                    if(dR34 < mindR34_debug && dR34 > 0.1)
                    {
                        numClosestPass34++;
                    }
                    if(dR34 < mindR34 && dR34 > 0.1)
                    {
                        mindR34 = dR34;
                        closestJet34 = &j;
                    }
                }
                //std::cout <<" mindR12 = "<< mindR12 << std::endl;
                //std::cout <<" mindR34 = "<< mindR34 << std::endl;
                //std::cout <<" mindR12_to1 = "<< mindR12_to1 << std::endl;
                //std::cout <<" mindR34_to1 = "<< mindR34_to1 << std::endl;
                //std::cout <<" mindR12 = "<< mindR12_debug << std::endl;
                //std::cout <<" mindR34 = "<< mindR34_debug << std::endl;
                //std::cout <<" numClosestPass12 = "<< numClosestPass12 << std::endl;
                //std::cout <<" numClosestPass34 = "<< numClosestPass34 << std::endl;
                _numClosestPass12 = numClosestPass12;
                _numClosestPass34 = numClosestPass34;
                if(closestJet12)
                {
                    _mW12 = (closestJet12->momentum() + leastTagged12.momentum()).mass();
                    _mt12 = (closestJet12->momentum() + b11.momentum() + b12.momentum()).mass();	
                    _dRW12 = mindR12;
                    closestJets.push_back(TagJet(*closestJet12,-99,-100));
                    _pt3rdJet12 = closestJet12->pt();
                    _eta3rdJet12 = closestJet12->eta();
                    _phi3rdJet12 = closestJet12->phi();
                    _y3rdJet12 = closestJet12->rap();
                    _m3rdJet12 = closestJet12->mass();
                }
                if(closestJet34)
                {
                    _mW34 = (closestJet34->momentum() + leastTagged34.momentum()).mass();
                    _mt34 = (closestJet34->momentum() + b21.momentum() + b22.momentum()).mass();	
                    _dRW34 = mindR34;
                    closestJets.push_back(TagJet(*closestJet34,-99,-101));
                    _pt3rdJet34 = closestJet34->pt();
                    _eta3rdJet34 = closestJet34->eta();
                    _phi3rdJet34 = closestJet34->phi();
                    _y3rdJet34 = closestJet34->rap();
                    _m3rdJet34 = closestJet34->mass();
                }
                if (closestJet12 || closestJet34) {
                    _sameClosest = 2;
                    //std::cout << "at least 1 closest Jet found!" << std::endl;
                    
                }
                if (closestJet12 && closestJet34) {
                    _sameClosest = 0;
                    if(closestJet12->eta() == closestJet34->eta() && closestJet12->phi() == closestJet34->phi()){
                        //std::cout << "closest Jets are same jets!" << std::endl;
                        _sameClosest = 1;
                        //std::cout <<" pt = "<< closestJet12->pt() <<" GeV, eta = "<< closestJet12->eta() <<", phi = "<< closestJet12->phi() << std::endl;
                        //std::cout <<" pt = "<< closestJet34->pt() <<" GeV, eta = "<< closestJet34->eta() <<", phi = "<< closestJet34->phi() << std::endl;
                        return;
                    }
                }
                if (_quarkMatchStatus == 3 && _alignStatus == 2) {
                    //std::cout << "This event is golden!" << std::endl;
                    double multiMatch3rd = 0;
                    double num3rdJetsMatched = 0;
                    _numClosestFound = closestJets.size();
                    std::map<int, quark_with_dR> quark;
                    if (closestJets.size() > 0) {
                        //std::cout << "some closest Jets Found!" << std::endl;
                        ParticleVector thirdPar;
                        thirdPar.push_back(tDaughters[thirdIndex_A]);
                        thirdPar.push_back(tDaughters[thirdIndex_B]);
                        int numToMatch = closestJets.size();
                        num3rdJetsMatched = matchJetToPar(numToMatch, multiMatch3rd, quark, closestJets, thirdPar);
                        _multiMatch3rd = multiMatch3rd;
                        _num3rdJetsMatched = num3rdJetsMatched;
                    }
                    if (closestJet12) {
                        _dR3rdJet12 = deltaR(*closestJet12, tDaughters[thirdIndex_A]);
                        //std::cout << "_dR3rdJet12: " << _dR3rdJet12 << std::endl;
                        map<int, quark_with_dR>::iterator it;
                        for(it = quark.begin(); it != quark.end(); it++) {
                            //std::cout << "key: " << it.first << ", dR: " << it.second.dR << std::endl;
                            if (it->second.dR == _dR3rdJet12) {
                                _quarkMatched3rdJet12 = it->second.quark;
                                _dRMatched3rdJet12 = it->second.dR;
                            }
                        }
                    }
                    if (closestJet34) {
                        _dR3rdJet34 = deltaR(*closestJet34, tDaughters[thirdIndex_B]);
                        //std::cout << "_dR3rdJet34: " << _dR3rdJet34 << std::endl;
                        map<int, quark_with_dR>::iterator it;
                        for(it = quark.begin(); it != quark.end(); it++) {
                            //std::cout << "key: " << it.first << ", dR: " << it.second.dR << std::endl;
                            if (it->second.dR == _dR3rdJet34) {
                                _quarkMatched3rdJet34 = it->second.quark;
                                _dRMatched3rdJet34 = it->second.dR;
                            }
                        }
                    }
                    if (num3rdJetsMatched == -5) {
                        _3rdMatchStatus = -5;
                    }
                    else if (closestJets.size() == 2 && num3rdJetsMatched <= 2) {
                        if (closestJets.size() == 2 && num3rdJetsMatched == 0) {
                            _3rdMatchStatus = -1;
                        }
                        else if (closestJets.size() == 2 && num3rdJetsMatched == 2) {
                            _3rdMatchStatus = 3;
                        }
                        else {
                            _3rdMatchStatus = 2;
                        }
                    }
                    else if (closestJets.size() == 1 && num3rdJetsMatched < 2) {
                        if (closestJets.size() == 1 && num3rdJetsMatched == 0) {
                            _3rdMatchStatus = -2;
                        }
                        else {
                            _3rdMatchStatus = 1;
                        }
                    }
                    else if (closestJets.size() == 0 && num3rdJetsMatched == 0) {
                        _3rdMatchStatus = 0;
                    }
                    else {
                        _3rdMatchStatus = -3;
                    }
                    /*
                    foreach(const TagJet& j, closestJets)
                    {
                        //std::cout << "RUNNING FILLER " << std::endl;
                        if (j.flavour() == -100) {
                            _pt3rdJet12 = j.pt();
                            _eta3rdJet12 = j.eta();
                            _phi3rdJet12 = j.phi();
                            _y3rdJet12 = j.rap();
                            _m3rdJet12 = j.mass();
                        }
                        if (j.flavour() == -101) {
                            _pt3rdJet34 = j.pt();
                            _eta3rdJet34 = j.eta();
                            _phi3rdJet34 = j.phi();
                            _y3rdJet34 = j.rap();
                            _m3rdJet34 = j.mass();
                        }
                    }
                    */
                    //std::cout << "3rdMatchStatus: " << _3rdMatchStatus << std::endl;

                }
            }
            /*
            /// Fill histograms for the 5 unique angles in H->bb [Phys. Rev. D 86 095031 (2012)]
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
             */
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
            /*
            const HepMC::GenParticle *findLastInChain(const HepMC::GenParticle *particle) {
                HepMC::GenVertex *decay = particle->end_vertex();
                std::cout << "findLastInChain Running" << std::endl;
                bool continueChain = false;
                const HepMC::GenParticle *nextInChain;
                if (decay) {
                    std::cout << "Number of decay products: " << decay->particles_out_size() << std::endl;
                    HepMC::GenVertex::particles_out_const_iterator it;
                    for (it=decay->particles_out_const_begin(); it!=decay->particles_out_const_end(); ++it) {
                        Particle p(*it);
                        std::cout << "  ID: " << p.pid() << " mass: " << p.mass() << " GeV, pt: " << p.pt() << " GeV, eta: " << p.eta() << ", phi: " << p.phi() << std::endl;
                        //std::cout << "  Momentum: " << p.p3() << std::endl;
                        if ((*it)->pdg_id() == particle->pdg_id() && (*it)->barcode() > particle->barcode()) {
                            continueChain = true;
                            nextInChain = *it;
                        }
                    }
                }
                if (continueChain) {
                    return findLastInChain(nextInChain);
                }
                else {
                    std::cout << "Thats the end of chain!" << std::endl;
                    return particle;
                }
            }
            */
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
				boost::shared_ptr<YODA::Histo1D> hist_pass = _histograms_1d[label+pass];
				boost::shared_ptr<YODA::Histo1D> hist_fail = _histograms_1d[label+fail];

				if (!hist_pass) return;

				for (int i=0; i<hist_pass->bins().size(); ++i) {
					//double entries_pass = hist_pass->bin(i).sumW();
					//double entries_total = entries_pass + (hist_fail ? hist_fail->bin(i).sumW() : 0.);

					/*if (entries_total > 0 && entries_pass >= 0) {
					  bookAndFill1D(label+"efficiency", hist_pass->title(), (hist_pass->bin(i).lowEdge()+hist_pass->bin(i).highEdge())/2., 
					  axis, std::string("BLARG"), hist_pass->lowEdge(), hist_pass->highEdge(), 

					//            axis, hist_pass->bins()->title(), hist_pass->lowEdge(), hist_pass->highEdge(), 
					entries_pass/entries_total * hist_pass->bin(i).width());
					}*/
				}
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

			bool dijetSelection(Jets &bjets) 
			{

				Jets input(bjets);
				bjets.clear();

				for (unsigned int i=0; i<input.size(); ++i) {
					for (unsigned int j=i+1; j<input.size(); ++j) {
						if (deltaR(input[i], input[j]) > 1.5) continue;
						FourMomentum dijet = input[i].momentum() + input[j].momentum();

						if (dijet.pT() < 150*GeV) continue;
						if (fabs(dijet.mass()/GeV - 115) > 30) continue;

						bjets.push_back(input[i]);
						bjets.push_back(input[j]);

						input.erase(input.begin() + j--);
						input.erase(input.begin() + i--);
						break;
					}
				}

				return bjets.size() >= 4;
			}
        	static bool tagGreater(const TagJet& a, const TagJet& b){ return (a.tagEff() > b.tagEff()); }

			bool dijetSelection(std::vector<TagJet> &bjets) 
			{
				std::vector<TagJet> input = bjets;	
				std::stable_sort(input.begin(), input.end(), tagGreater);
				/*for(int j = 0; j < input.size(); ++j)
				{
					std::cout<<"dijetSelection: "<< j <<", tagEff = "<< input[j].tagEff() <<", pt = "<< input[j].pt() << std::endl;
				}*/
				bjets.clear();

				/*for (unsigned int i=0; i<4; ++i) 
				{
					for (unsigned int j=i+1; j<4; ++j) 
					{*/
				for (unsigned int i=0; i<input.size(); ++i) 
				{
					for (unsigned int j=i+1; j<input.size(); ++j) 
					{
						if (deltaR(input[i], input[j]) > 1.5) continue;
						FourMomentum dijet = input[i].momentum() + input[j].momentum();

						if (dijet.pT() < 150*GeV) continue;
						//if (fabs(dijet.mass()/GeV - 115) > 25) continue;

						bjets.push_back(input[i]);
						bjets.push_back(input[j]);

						input.erase(input.begin() + j--);
						input.erase(input.begin() + i--);
						break;
					}
				}

				return bjets.size() >= 4;
			}

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
                
                //std::cout << (q1+q2).pT() << std::endl;

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
        
            double _ptA, _etaA, _mA, _yA, _phiA;
            double _ptB, _etaB, _mB, _yB, _phiB;
        
            double _quarkMiss12, _ptMiss12, _etaMiss12, _phiMiss12, _yMiss12, _mMiss12;
            double _quarkMiss34, _ptMiss34, _etaMiss34, _phiMiss34, _yMiss34, _mMiss34;
            double _numJetsMatched, _numQuarks, _quarkMatchStatus, _alignStatus, _multiMatch;
            double _numNonDup, _ptMinAll, _etaMaxAll, _numClosestFound, _multiMatch3rd, _num3rdJetsMatched, _sameClosest, _3rdMatchStatus;
            double _charge12, _charge34;
            double _bLeadStatus12, _bLeadStatus34;
            double _dijetsHaveTwoBquarks, _bothBleading;
            double _dR_sublead12, _dR_sublead34;
            double _dR_lead12, _dR_lead34;
            double _dR_Wquarks12, _dR_Wquarks34;
            double _quark1, _quark2, _quark3, _quark4;
            double _deltaR1, _deltaR2, _deltaR3, _deltaR4;
            double thirdIndex_A, thirdIndex_B;
            double _numClosestPass12, _numClosestPass34;
        
            double _quarkMatched3rdJet12, _pt3rdJet12, _eta3rdJet12, _phi3rdJet12, _y3rdJet12, _m3rdJet12, _dR3rdJet12, _dRMatched3rdJet12;
            double _quarkMatched3rdJet34, _pt3rdJet34, _eta3rdJet34, _phi3rdJet34, _y3rdJet34, _m3rdJet34, _dR3rdJet34, _dRMatched3rdJet34;

			/// @name Histograms
			//@{
			map<string, shared_ptr<YODA::Histo1D> > _histograms_1d;
			map<string, shared_ptr<YODA::Histo2D> > _histograms_2d;
			//@}


	};



	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(Summer_2014_Study);

}
