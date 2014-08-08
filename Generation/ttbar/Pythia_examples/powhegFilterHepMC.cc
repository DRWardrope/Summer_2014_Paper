// main31.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Richard Corke, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

//#include "FastJetFilter.h"
#include "Pythia8/FastJet3.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
using namespace Pythia8;

//==========================================================================

// Use userhooks to veto PYTHIA emissions above the POWHEG scale.
class PowhegFilterHooks : public UserHooks {

	public:

		// Constructor and destructor.
		PowhegFilterHooks(int nFinalIn, int vetoModeIn, int vetoCountIn,
				int pThardModeIn, int pTemtModeIn, int emittedModeIn,
				int pTdefModeIn, int MPIvetoModeIn) : nFinal(nFinalIn),
		vetoMode(vetoModeIn), vetoCount(vetoCountIn),
		pThardMode(pThardModeIn), pTemtMode(pTemtModeIn),
		emittedMode(emittedModeIn), pTdefMode(pTdefModeIn),
		MPIvetoMode(MPIvetoModeIn) {};
		~PowhegFilterHooks() {}

		//Filtering hooks
		bool canVetoProcessLevel() { return false;}
		bool doVetoProcessLevel(Event& process)
		{
			//std::cout<<"Hello from doVetoProcessLevel, there are "<< process.size() <<" particles in this event!"<< std::endl;
			vector<float> tagProb;
			for(int i = 0; i < process.size(); ++i)
			{
				if(process[i].status() != 23) continue;
				if(process[i].pT() < 30.) continue;
				if(fabs(process[i].eta()) > 2.7) continue;
				if(abs(process[i].id()) == 16) continue;
				//taggable.push_back(abs(process[i]->id()));
				//std::cout<<"Particle "<< i <<", status = "<< process[i].status() <<", pdgID = "<< process[i].id() <<", mother1 = "<< process[i].mother1() <<", p = "<< process[i].p() << std::endl;
				//std::cout<<"Particle "<< i <<", status = "<< process[i].status() <<", pdgID = "<< process[i].id() <<", pt = "<< process[i].pT() <<"GeV, eta = "<< process[i].eta() << std::endl;
				int id1 = abs(process[i].id());
				if(id1 == 5) tagProb.push_back(0.7);
				else if(id1 == 4 || id1 == 15) tagProb.push_back(0.2);
				else if(id1 == 21)
				{
					tagProb.push_back(0.0135); //P(g->bb)*(1-eff(b)^2) + P(g->cc)*(1 - eff(c)^2)
					//tagProb.push_back(0.008039); //P(g->bb)*eff(b) + P(g->cc)*eff(c)
					//tagProb.push_back(0.008039); //Based on measurements in PDG, end of chap. 19.
				}else tagProb.push_back(0.01);
			}
			if(tagProb.size() < 4) return true;
			//else return false;

			float evtTagProb = 0.;
			for(std::vector<float>::const_iterator p1 = tagProb.begin(); p1 != tagProb.end(); ++p1)
			{
				for(std::vector<float>::const_iterator p2 = tagProb.begin(); p2 != tagProb.end(); ++p2)
				{
					if(p2 == p1) continue;
					for(std::vector<float>::const_iterator p3 = tagProb.begin(); p3 != tagProb.end(); ++p3)
					{
						if(p3 == p2 || p3 == p1) continue;
						for(std::vector<float>::const_iterator p4 = tagProb.begin(); p4 != tagProb.end(); ++p4)
						{
							if(p4 == p3 || p4 == p2 || p4 == p1) continue;
							evtTagProb += (*p1)*(*p2)*(*p3)*(*p4); 
						}
					}
				}
			}
			//std::cout<<"Event tag probability = "<< evtTagProb << std::endl;
			if(evtTagProb > 0.05) return false;
			else return true;
		}
		//--------------------------------------------------------------------------

		// Routines to calculate the pT (according to pTdefMode) in a splitting:
		//   ISR: i (radiator after)  -> j (emitted after) k (radiator before)
		//   FSR: i (radiator before) -> j (emitted after) k (radiator after)
		// For the Pythia pT definition, a recoiler (after) must be specified.

		// Compute the Pythia pT separation. Based on pTLund function in History.cc
		double pTpythia(const Event &e, int RadAfterBranch, int EmtAfterBranch,
				int RecAfterBranch, bool FSR) {

			// Convenient shorthands for later
			Vec4 radVec = e[RadAfterBranch].p();
			Vec4 emtVec = e[EmtAfterBranch].p();
			Vec4 recVec = e[RecAfterBranch].p();
			int  radID  = e[RadAfterBranch].id();

			// Calculate virtuality of splitting
			double sign = (FSR) ? 1. : -1.;
			Vec4 Q(radVec + sign * emtVec);
			double Qsq = sign * Q.m2Calc();

			// Mass term of radiator
			double m2Rad = (abs(radID) >= 4 && abs(radID) < 7) ?
				pow2(particleDataPtr->m0(radID)) : 0.;

			// z values for FSR and ISR
			double z, pTnow;
			if (FSR) {
				// Construct 2 -> 3 variables
				Vec4 sum = radVec + recVec + emtVec;
				double m2Dip = sum.m2Calc();
				double x1 = 2. * (sum * radVec) / m2Dip;
				double x3 = 2. * (sum * emtVec) / m2Dip;
				z     = x1 / (x1 + x3);
				pTnow = z * (1. - z);

			} else {
				// Construct dipoles before/after splitting
				Vec4 qBR(radVec - emtVec + recVec);
				Vec4 qAR(radVec + recVec);
				z     = qBR.m2Calc() / qAR.m2Calc();
				pTnow = (1. - z);
			}

			// Virtuality with correct sign
			pTnow *= (Qsq - sign * m2Rad);

			// Can get negative pT for massive splittings
			if (pTnow < 0.) {
				cout << "Warning: pTpythia was negative" << endl;
				return -1.;
			}

#ifdef DBGOUTPUT
			cout << "pTpythia: rad = " << RadAfterBranch << ", emt = "
				<< EmtAfterBranch << ", rec = " << RecAfterBranch
				<< ", pTnow = " << sqrt(pTnow) << endl;
#endif

			// Return pT
			return sqrt(pTnow);
		}

		// Compute the POWHEG pT separation between i and j
		double pTpowheg(const Event &e, int i, int j, bool FSR) {

			// pT value for FSR and ISR
			double pTnow = 0.;
			if (FSR) {
				// POWHEG d_ij (in CM frame). Note that the incoming beams have not
				// been updated in the parton systems pointer yet (i.e. prior to any
				// potential recoil).
				int iInA = partonSystemsPtr->getInA(0);
				int iInB = partonSystemsPtr->getInB(0);
				double betaZ = - ( e[iInA].pz() + e[iInB].pz() ) /
					( e[iInA].e()  + e[iInB].e()  );
				Vec4 iVecBst(e[i].p()), jVecBst(e[j].p());
				iVecBst.bst(0., 0., betaZ);
				jVecBst.bst(0., 0., betaZ);
				pTnow = sqrt( (iVecBst + jVecBst).m2Calc() *
						iVecBst.e() * jVecBst.e() /
						pow2(iVecBst.e() + jVecBst.e()) );

			} else {
				// POWHEG pT_ISR is just kinematic pT
				pTnow = e[j].pT();
			}

			// Check result
			if (pTnow < 0.) {
				cout << "Warning: pTpowheg was negative" << endl;
				return -1.;
			}

#ifdef DBGOUTPUT
			cout << "pTpowheg: i = " << i << ", j = " << j
				<< ", pTnow = " << pTnow << endl;
#endif

			return pTnow;
		}

		// Calculate pT for a splitting based on pTdefMode.
		// If j is -1, all final-state partons are tried.
		// If i, k, r and xSR are -1, then all incoming and outgoing
		// partons are tried.
		// xSR set to 0 means ISR, while xSR set to 1 means FSR
		double pTcalc(const Event &e, int i, int j, int k, int r, int xSRin) {

			// Loop over ISR and FSR if necessary
			double pTemt = -1., pTnow;
			int xSR1 = (xSRin == -1) ? 0 : xSRin;
			int xSR2 = (xSRin == -1) ? 2 : xSRin + 1;
			for (int xSR = xSR1; xSR < xSR2; xSR++) {
				// FSR flag
				bool FSR = (xSR == 0) ? false : true;

				// If all necessary arguments have been given, then directly calculate.
				// POWHEG ISR and FSR, need i and j.
				if ((pTdefMode == 0 || pTdefMode == 1) && i > 0 && j > 0) {
					pTemt = pTpowheg(e, i, j, (pTdefMode == 0) ? false : FSR);

					// Pythia ISR, need i, j and r.
				} else if (!FSR && pTdefMode == 2 && i > 0 && j > 0 && r > 0) {
					pTemt = pTpythia(e, i, j, r, FSR);

					// Pythia FSR, need k, j and r.
				} else if (FSR && pTdefMode == 2 && j > 0 && k > 0 && r > 0) {
					pTemt = pTpythia(e, k, j, r, FSR);

					// Otherwise need to try all possible combintations.
				} else {
					// Start by finding incoming legs to the hard system after
					// branching (radiator after branching, i for ISR).
					// Use partonSystemsPtr to find incoming just prior to the
					// branching and track mothers.
					int iInA = partonSystemsPtr->getInA(0);
					int iInB = partonSystemsPtr->getInB(0);
					while (e[iInA].mother1() != 1) { iInA = e[iInA].mother1(); }
					while (e[iInB].mother1() != 2) { iInB = e[iInB].mother1(); }

					// If we do not have j, then try all final-state partons
					int jNow = (j > 0) ? j : 0;
					int jMax = (j > 0) ? j + 1 : e.size();
					for (; jNow < jMax; jNow++) {

						// Final-state and coloured jNow only
						if (!e[jNow].isFinal() || e[jNow].colType() == 0) continue;

						// POWHEG
						if (pTdefMode == 0 || pTdefMode == 1) {

							// ISR - only done once as just kinematical pT
							if (!FSR) {
								pTnow = pTpowheg(e, iInA, jNow, (pTdefMode == 0) ? false : FSR);
								if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);

								// FSR - try all outgoing partons from system before branching
								// as i. Note that for the hard system, there is no
								// "before branching" information.
							} else {

								int outSize = partonSystemsPtr->sizeOut(0);
								for (int iMem = 0; iMem < outSize; iMem++) {
									int iNow = partonSystemsPtr->getOut(0, iMem);

									// Coloured only, i != jNow and no carbon copies
									if (iNow == jNow || e[iNow].colType() == 0) continue;
									if (jNow == e[iNow].daughter1()
											&& jNow == e[iNow].daughter2()) continue;

									pTnow = pTpowheg(e, iNow, jNow, (pTdefMode == 0)
											? false : FSR);
									if (pTnow > 0.) pTemt = (pTemt < 0)
										? pTnow : min(pTemt, pTnow);
								} // for (iMem)

							} // if (!FSR)

							// Pythia
						} else if (pTdefMode == 2) {

							// ISR - other incoming as recoiler
							if (!FSR) {
								pTnow = pTpythia(e, iInA, jNow, iInB, FSR);
								if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);
								pTnow = pTpythia(e, iInB, jNow, iInA, FSR);
								if (pTnow > 0.) pTemt = (pTemt < 0) ? pTnow : min(pTemt, pTnow);

								// FSR - try all final-state coloured partons as radiator
								//       after emission (k).
							} else {
								for (int kNow = 0; kNow < e.size(); kNow++) {
									if (kNow == jNow || !e[kNow].isFinal() ||
											e[kNow].colType() == 0) continue;

									// For this kNow, need to have a recoiler.
									// Try two incoming.
									pTnow = pTpythia(e, kNow, jNow, iInA, FSR);
									if (pTnow > 0.) pTemt = (pTemt < 0)
										? pTnow : min(pTemt, pTnow);
									pTnow = pTpythia(e, kNow, jNow, iInB, FSR);
									if (pTnow > 0.) pTemt = (pTemt < 0)
										? pTnow : min(pTemt, pTnow);

									// Try all other outgoing.
									for (int rNow = 0; rNow < e.size(); rNow++) {
										if (rNow == kNow || rNow == jNow ||
												!e[rNow].isFinal() || e[rNow].colType() == 0) continue;
										pTnow = pTpythia(e, kNow, jNow, rNow, FSR);
										if (pTnow > 0.) pTemt = (pTemt < 0)
											? pTnow : min(pTemt, pTnow);
									} // for (rNow)

								} // for (kNow)
							} // if (!FSR)
						} // if (pTdefMode)
					} // for (j)
				}
			} // for (xSR)

#ifdef DBGOUTPUT
			cout << "pTcalc: i = " << i << ", j = " << j << ", k = " << k
				<< ", r = " << r << ", xSR = " << xSRin
				<< ", pTemt = " << pTemt << endl;
#endif

			return pTemt;
		}

		//--------------------------------------------------------------------------

		// Extraction of pThard based on the incoming event.
		// Assume that all the final-state particles are in a continuous block
		// at the end of the event and the final entry is the POWHEG emission.
		// If there is no POWHEG emission, then pThard is set to SCALUP.

		bool canVetoMPIStep()    { return true; }
		int  numberVetoMPIStep() { return 1; }
		bool doVetoMPIStep(int nMPI, const Event &e) {
			// Extra check on nMPI
			if (nMPI > 1) return false;

			// Find if there is a POWHEG emission. Go backwards through the
			// event record until there is a non-final particle. Also sum pT and
			// find pT_1 for possible MPI vetoing
			int    count = 0;
			double pT1 = 0., pTsum = 0.;
			for (int i = e.size() - 1; i > 0; i--) {
				if (e[i].isFinal()) {
					count++;
					pT1    = e[i].pT();
					pTsum += e[i].pT();
				} else break;
			}
			// Extra check that we have the correct final state
			if (count != nFinal && count != nFinal + 1) {
				cout << "Error: wrong number of final state particles in event" << endl;
				exit(1);
			}
			// Flag if POWHEG radiation present and index
			bool isEmt = (count == nFinal) ? false : true;
			int  iEmt  = (isEmt) ? e.size() - 1 : -1;

			// If there is no radiation or if pThardMode is 0 then set pThard = SCALUP.
			if (!isEmt || pThardMode == 0) {
				pThard = infoPtr->scalup();

				// If pThardMode is 1 then the pT of the POWHEG emission is checked against
				// all other incoming and outgoing partons, with the minimal value taken
			} else if (pThardMode == 1) {
				pThard = pTcalc(e, -1, iEmt, -1, -1, -1);

				// If pThardMode is 2, then the pT of all final-state partons is checked
				// against all other incoming and outgoing partons, with the minimal value
				// taken
			} else if (pThardMode == 2) {
				pThard = pTcalc(e, -1, -1, -1, -1, -1);
			}

			// Find MPI veto pT if necessary
			if (MPIvetoMode == 1) {
				pTMPI = (isEmt) ? pTsum / 2. : pT1;
			}

#ifdef DBGOUTPUT
			cout << "doVetoMPIStep: Qfac = " << infoPtr->scalup()
				<< ", pThard = " << pThard << endl << endl;
#endif

			// Initialise other variables
			accepted   = false;
			nAcceptSeq = nISRveto = nFSRveto = 0;

			// Do not veto the event
			return false;
		}

		//--------------------------------------------------------------------------

		// ISR veto

		bool canVetoISREmission() { return (vetoMode == 0) ? false : true; }
		bool doVetoISREmission(int, const Event &e, int iSys) {
			// Must be radiation from the hard system
			if (iSys != 0) return false;

			// If we already have accepted 'vetoCount' emissions in a row, do nothing
			if (vetoMode == 1 && nAcceptSeq >= vetoCount) return false;

			// Pythia radiator after, emitted and recoiler after.
			int iRadAft = -1, iEmt = -1, iRecAft = -1;
			for (int i = e.size() - 1; i > 0; i--) {
				if      (iRadAft == -1 && e[i].status() == -41) iRadAft = i;
				else if (iEmt    == -1 && e[i].status() ==  43) iEmt    = i;
				else if (iRecAft == -1 && e[i].status() == -42) iRecAft = i;
				if (iRadAft != -1 && iEmt != -1 && iRecAft != -1) break;
			}
			if (iRadAft == -1 || iEmt == -1 || iRecAft == -1) {
				e.list();
				cout << "Error: couldn't find Pythia ISR emission" << endl;
				exit(1);
			}

			// pTemtMode == 0: pT of emitted w.r.t. radiator
			// pTemtMode == 1: min(pT of emitted w.r.t. all incoming/outgoing)
			// pTemtMode == 2: min(pT of all outgoing w.r.t. all incoming/outgoing)
			int xSR      = (pTemtMode == 0) ? 0       : -1;
			int i        = (pTemtMode == 0) ? iRadAft : -1;
			int j        = (pTemtMode != 2) ? iEmt    : -1;
			int k        = -1;
			int r        = (pTemtMode == 0) ? iRecAft : -1;
			double pTemt = pTcalc(e, i, j, k, r, xSR);

#ifdef DBGOUTPUT
			cout << "doVetoISREmission: pTemt = " << pTemt << endl << endl;
#endif

			// Veto if pTemt > pThard
			if (pTemt > pThard) {
				nAcceptSeq = 0;
				nISRveto++;
				return true;
			}

			// Else mark that an emission has been accepted and continue
			nAcceptSeq++;
			accepted = true;
			return false;
		}

		//--------------------------------------------------------------------------

		// FSR veto

		bool canVetoFSREmission() { return (vetoMode == 0) ? false : true; }
		bool doVetoFSREmission(int, const Event &e, int iSys, bool) {
			// Must be radiation from the hard system
			if (iSys != 0) return false;

			// If we already have accepted 'vetoCount' emissions in a row, do nothing
			if (vetoMode == 1 && nAcceptSeq >= vetoCount) return false;

			// Pythia radiator (before and after), emitted and recoiler (after)
			int iRecAft = e.size() - 1;
			int iEmt    = e.size() - 2;
			int iRadAft = e.size() - 3;
			int iRadBef = e[iEmt].mother1();
			if ( (e[iRecAft].status() != 52 && e[iRecAft].status() != -53) ||
					e[iEmt].status() != 51 || e[iRadAft].status() != 51) {
				e.list();
				cout << "Error: couldn't find Pythia FSR emission" << endl;
				exit(1);
			}

			// Behaviour based on pTemtMode:
			//  0 - pT of emitted w.r.t. radiator before
			//  1 - min(pT of emitted w.r.t. all incoming/outgoing)
			//  2 - min(pT of all outgoing w.r.t. all incoming/outgoing)
			int xSR = (pTemtMode == 0) ? 1       : -1;
			int i   = (pTemtMode == 0) ? iRadBef : -1;
			int k   = (pTemtMode == 0) ? iRadAft : -1;
			int r   = (pTemtMode == 0) ? iRecAft : -1;

			// When pTemtMode is 0 or 1, iEmt has been selected
			double pTemt = 0.;
			if (pTemtMode == 0 || pTemtMode == 1) {
				// Which parton is emitted, based on emittedMode:
				//  0 - Pythia definition of emitted
				//  1 - Pythia definition of radiated after emission
				//  2 - Random selection of emitted or radiated after emission
				//  3 - Try both emitted and radiated after emission
				int j = iRadAft;
				if (emittedMode == 0 || (emittedMode == 2 && rndmPtr->flat() < 0.5)) j++;

				for (int jLoop = 0; jLoop < 2; jLoop++) {
					if      (jLoop == 0) pTemt = pTcalc(e, i, j, k, r, xSR);
					else if (jLoop == 1) pTemt = min(pTemt, pTcalc(e, i, j, k, r, xSR));

					// For emittedMode == 3, have tried iRadAft, now try iEmt
					if (emittedMode != 3) break;
					if (k != -1) swap(j, k); else j = iEmt;
				}

				// If pTemtMode is 2, then try all final-state partons as emitted
			} else if (pTemtMode == 2) {
				pTemt = pTcalc(e, i, -1, k, r, xSR);

			}

#ifdef DBGOUTPUT
			cout << "doVetoFSREmission: pTemt = " << pTemt << endl << endl;
#endif

			// Veto if pTemt > pThard
			if (pTemt > pThard) {
				nAcceptSeq = 0;
				nFSRveto++;
				return true;
			}

			// Else mark that an emission has been accepted and continue
			nAcceptSeq++;
			accepted = true;
			return false;
		}

		//--------------------------------------------------------------------------

		// MPI veto

		bool canVetoMPIEmission() { return (MPIvetoMode == 0) ? false : true; }
		bool doVetoMPIEmission(int, const Event &e) {
			if (MPIvetoMode == 1) {
				if (e[e.size() - 1].pT() > pTMPI) {
#ifdef DBGOUTPUT
					cout << "doVetoMPIEmission: pTnow = " << e[e.size() - 1].pT()
						<< ", pTMPI = " << pTMPI << endl << endl;
#endif
					return true;
				}
			}
			return false;
		}

		//--------------------------------------------------------------------------

		// Functions to return information

		int    getNISRveto() { return nISRveto; }
		int    getNFSRveto() { return nFSRveto; }

	private:
		int    nFinal, vetoMode, vetoCount, pThardMode, pTemtMode,
			   emittedMode, pTdefMode, MPIvetoMode;
		double pThard, pTMPI;
		bool   accepted;
		// The number of accepted emissions (in a row)
		int nAcceptSeq;
		// Statistics on vetos
		unsigned long int nISRveto, nFSRveto;

};

class FastJetFilter{
	public:
		bool isStableBHadron(int pdg) 
		{
			pdg=abs(pdg);
			// For multiplicity only consider ground states
			return ((pdg==511 || pdg==521 || pdg==531 || pdg==541 || pdg==5122) ||
					(pdg==5132 || pdg==5142 || pdg==5232 || pdg==5242) ||
					(pdg==5332 || pdg==5342 || pdg==5432 || pdg==5434) ||
					(pdg==5442 || pdg==5444 || pdg==5512 || pdg==5514) ||
					(pdg==5522 || pdg==5524 || pdg==5532 || pdg==5542) ||
					(pdg==5544 || pdg==5554 ));
		}
		bool isCHadron(int pdg) 
		{
			pdg=abs(pdg);  
			return ((pdg>=400&&pdg<500&&pdg !=443)|| (pdg>=4000 && pdg<5000)
					|| (pdg>=10411 && pdg<=10445) || (pdg>=20411 && pdg<=20445));
		}

		bool reject(Event& event)
		{
			// Fastjet analysis - select algorithm and parameters
			fastjet::Strategy              strategy = fastjet::Best;
			fastjet::RecombinationScheme   recombScheme = fastjet::E_scheme;
			fastjet::JetDefinition         jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4,recombScheme, strategy);

			// Make Fastjet particles
			std::vector<fastjet::PseudoJet> parts;
			std::vector<fastjet::PseudoJet> taggables;
			fastjet::Selector visSel = fastjet::SelectorIsVisible();
			fastjet::Selector finalSel = fastjet::SelectorIsFinal();
			for(int i = 0; i < event.size(); ++i)
			{
				if(fabs(event[i].eta()) > 2.9) continue;
				fastjet::PseudoJet fj_particle = event[i];
				if(visSel.pass(fj_particle) && finalSel.pass(fj_particle)) parts.push_back(fj_particle);	
			/*	if(isStableBHadron(event[i].id()) || isCHadron(event[i].id()) || abs(event[i].id()) == 15) 
				{
					//taggables.push_back((0.001/fj_particle.modp())*fj_particle);
					parts.push_back((0.001/fj_particle.modp())*fj_particle);
				}*/
			}
			//Select particles for jet-clustering
//			std::cout<<"FastJetFilter::reject: There are "<< parts.size() <<" visible, final state particles within |eta| < 2.9"<< std::endl;
			/*std::cout<<"FastJetFilter::reject: There are "<< taggables.size() <<" taggable hadrons within |eta| < 2.9"<< std::endl;
			for(std::vector<fastjet::PseudoJet>::const_iterator taggable = taggables.begin(); taggable != taggables.end(); ++taggable)
			{
				std::cout<<"Taggable: pdg = "<< taggable->user_info<Particle>().id() <<", status = "<< taggable->user_info<Particle>().status() 
							<<", pt = "<< taggable->pt() <<" GeV, E = "<< taggable->E() <<" GeV"<< std::endl;
			}*/

			//Jet Clustering
			fastjet::ClusterSequence clustSeq(parts, jetDef);
			std::vector<fastjet::PseudoJet> looseJets = clustSeq.inclusive_jets(40.);//sorted_by_pt(clustSeq.inclusive_jets());
			fastjet::Selector absEtaSel = fastjet::SelectorAbsEtaMax(2.5);
			std::vector<fastjet::PseudoJet> jets = absEtaSel(looseJets);
//			std::cout<<"FastJetFilter::reject has found "<< jets.size() <<" jets."<< std::endl;
			//Decision-making	
			if(jets.size() < 4)
			{
//				std::cout<<"Reject the event!"<< std::endl;
				return true;
			}else return false;
		}
	private:
};
//==========================================================================
			/*}else{
				int nB = 0;
				for(std::vector<fastjet::PseudoJet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet)
				{
					std::cout<<"Found jet pt = "<< jet->pt() <<", eta = "<< jet->eta() << std::endl;
					std::vector<fastjet::PseudoJet> constituents = jet->constituents();
					int jetPDG = 0;
					for(std::vector<fastjet::PseudoJet>::const_iterator con = constituents.begin(); con != constituents.end(); ++con)
					{
					  	std::cout<<"----> PDG = "<< con->user_info<Particle>().id() <<", status = "<< con->user_info<Particle>().status() << std::endl;
						if(isStableBHadron(con->user_info<Particle>().id()))
						{
							jetPDG = 5;	
							++nB;
							break;
						}
						else if(isCHadron(con->user_info<Particle>().id()) || abs(con->user_info<Particle>().id()) == 15) //Treat taus as charms for sake of filtering
						{
							jetPDG = 4;	
						}
					}
					std::cout<<"Jet flavour = "<< jetPDG << std::endl;
				}
				std::cout<<"Keep the event!"<< std::endl;	
				if(nB == 0) return true;
				else return false;
			}*/

int main(int argc, char* argv[]) 
{
	// Check that correct number of command-line arguments
	if (argc != 3) {
		cerr << " Unexpected number of command-line arguments. \n You are"
			<< " expected to provide one input and one output file name. \n"
			<< " Program stopped! " << endl;
		return 1;
	}

	// Check that the provided input name corresponds to an existing file.
	ifstream is(argv[1]);
	if (!is) {
		cerr << " Command-line file " << argv[1] << " was not found. \n"
			<< " Program stopped! " << endl;
		return 1;
	}

	// Confirm that external files will be used for input and output.
	cout << "\n >>> PYTHIA settings will be read from file " << argv[1]
		<< " <<< \n >>> HepMC events will be written to file "
		<< argv[2] << " <<< \n" << endl;

	// Interface for conversion from Pythia8::Event to HepMC event.
	HepMC::Pythia8ToHepMC ToHepMC;

	// Specify file where HepMC events will be stored.
	HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

	// Generator
	Pythia pythia;
	// Add further settings that can be set in the configuration file
	pythia.settings.addMode("POWHEG:nFinal",    2, true, false, 1, 0);
	pythia.settings.addMode("POWHEG:veto",      0, true, true,  0, 2);
	pythia.settings.addMode("POWHEG:vetoCount", 3, true, false, 0, 0);
	pythia.settings.addMode("POWHEG:pThard",    0, true, true,  0, 2);
	pythia.settings.addMode("POWHEG:pTemt",     0, true, true,  0, 2);
	pythia.settings.addMode("POWHEG:emitted",   0, true, true,  0, 3);
	pythia.settings.addMode("POWHEG:pTdef",     0, true, true,  0, 2);
	pythia.settings.addMode("POWHEG:MPIveto",   0, true, true,  0, 1);
	// Load configuration file
	pythia.readFile( argv[1]);



	// Read in main settings
	int nEvent      = pythia.settings.mode("Main:numberOfEvents");
	int nError      = pythia.settings.mode("Main:timesAllowErrors");
	// Read in POWHEG settings
	int nFinal      = pythia.settings.mode("POWHEG:nFinal");
	int vetoMode    = pythia.settings.mode("POWHEG:veto");
	int vetoCount   = pythia.settings.mode("POWHEG:vetoCount");
	int pThardMode  = pythia.settings.mode("POWHEG:pThard");
	int pTemtMode   = pythia.settings.mode("POWHEG:pTemt");
	int emittedMode = pythia.settings.mode("POWHEG:emitted");
	int pTdefMode   = pythia.settings.mode("POWHEG:pTdef");
	int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
	bool loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);

	// Add in user hooks for shower vetoing
	PowhegFilterHooks *powhegHooks = NULL;
	if (loadHooks) {

		// Set ISR and FSR to start at the kinematical limit
		if (vetoMode > 0) {
			pythia.readString("SpaceShower:pTmaxMatch = 2");
			pythia.readString("TimeShower:pTmaxMatch = 2");
		}

		// Set MPI to start at the kinematical limit
		if (MPIvetoMode > 0) {
			pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
		}

		powhegHooks = new PowhegFilterHooks(nFinal, vetoMode, vetoCount,
				pThardMode, pTemtMode, emittedMode, pTdefMode, MPIvetoMode);
		pythia.setUserHooksPtr((UserHooks *) powhegHooks);
	}

	// Initialise and list settings
	pythia.init();
	//Initialize FastjetFilter object
	FastJetFilter fjFilter;

	// Counters for number of ISR/FSR emissions vetoed
	unsigned long int nISRveto = 0, nFSRveto = 0;

	// Begin event loop; generate until nEvent events are processed
	// or end of LHEF file
	int iEvent = 0, iError = 0;
	while (true) {

		// Generate the next event
		if (!pythia.next()) {

			// If failure because reached end of file then exit event loop
			if (pythia.info.atEndOfFile()) break;

			// Otherwise count event failure and continue/exit as necessary
			cout << "Warning: event " << iEvent << " failed" << endl;
			if (++iError == nError) {
				cout << "Error: too many event failures.. exiting" << endl;
				break;
			}
			continue;
		}

		/*
		 * Process dependent checks and analysis may be inserted here
		 */

		// Update ISR/FSR veto counters
		if (loadHooks) {
			nISRveto += powhegHooks->getNISRveto();
			nFSRveto += powhegHooks->getNFSRveto();
		}

		if(fjFilter.reject(pythia.event)) continue;
		// Construct new empty HepMC event and fill it.
		// Units will be as chosen for HepMC build, but can be changed
		// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
		HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
		ToHepMC.fill_next_event( pythia, hepmcevt );

		// Write the HepMC event to file. Done with it.
		//  ascii_io<<"Hello from event " << iEvent << std::endl;
		ascii_io << hepmcevt;
		delete hepmcevt;

		// If nEvent is set, check and exit loop if necessary
		++iEvent;
		if (nEvent != 0 && iEvent == nEvent) break;

	} // End of event loop.

	// Statistics, histograms and veto information
	pythia.stat();
	cout << "Number of ISR emissions vetoed: " << nISRveto << endl;
	cout << "Number of FSR emissions vetoed: " << nFSRveto << endl;
	cout << endl;

	// Done.
	if (powhegHooks) delete powhegHooks;
	return 0;
}

