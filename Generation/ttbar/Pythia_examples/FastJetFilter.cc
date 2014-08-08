// FastJetFilter.cc 
// Author: David Wardrope, David.Wardrope@cern.ch
// Decide whether event should be written out or not

#include "FastJetFilter.h"

bool FastJetFilter::reject(Event& event)
{
	// Fastjet analysis - select algorithm and parameters
	fastjet::Strategy              strategy = fastjet::Best;
	fastjet::RecombinationScheme   recombScheme = fastjet::E_scheme;
	fastjet::JetDefinition         jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4,recombScheme, strategy);

	// Make Fastjet particles
	std::vector<fastjet::PseudoJet> parts;
	for(int i = 0; i < event.size(); ++i)
	{
		//fastjet::PseudoJet fj_particle = event[i];
		parts.push_back(event[i]);	
	}
	//Select particles for jet-clustering
	
	//Jet Clustering
	fastjet::ClusterSequence clustSeq(parts, jetDef);
	std::vector<fastjet::PseudoJet> jets = sorted_by_pt(clustSeq.inclusive_jets());
	//Decision-making	
	if(jets.size() < 4) return true;
	else return false;
}
