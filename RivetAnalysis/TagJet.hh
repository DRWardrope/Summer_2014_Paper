#ifndef TAGJET_HH
#define TAGJET_HH
#include "Rivet/Jet.hh"


class TagJet : public Rivet::Jet
{
	public:
		TagJet(): Rivet::Jet(), _tagEff(-99.), _flav(-99) {};
		TagJet(Rivet::Jet& j, double tagEff=-99, int flav=-99): Rivet::Jet(j), _tagEff(tagEff), _flav(flav) {};
		TagJet(const Rivet::Jet& j, double tagEff=-99, int flav=-99): Rivet::Jet(j), _tagEff(tagEff), _flav(flav) {};

		int flavour() const { return _flav; }
		double tagEff() const { return _tagEff; }
		void setTagEff(double tagEff) { _tagEff = tagEff; }
		void setFlavour(int flav) { _flav = flav; }
	private:
		int _flav;
		double _tagEff;
};
#endif
