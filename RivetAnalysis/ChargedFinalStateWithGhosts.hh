#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

/**
 * ChargedFinalState projection that allows injecting "momentumless" particles into
 * the jet clustering. This way you can unambiguously associate hard process
 * quarks to particle jets. Ghost particles are specified by PDG ID.
 */
class ChargedFinalStateWithGhosts : public ChargedFinalState {

    public:
        
        /**
         * Constructor. Accepts an eta range and a minimal pT.
         */
        ChargedFinalStateWithGhosts(double mineta = -100000,
                double maxeta = 100000,
                double minpt = 0.0*GeV) 
            : ChargedFinalState(mineta, maxeta, minpt) {}

        /**
         * Clone this projection.
         */
        virtual const Projection *clone() const {
            return new ChargedFinalStateWithGhosts(*this);
        }

        /**
         * Add the specified pid and -pid to the list of ghost particles to
         * include.
         */
        ChargedFinalStateWithGhosts &ghostIdPair(PdgId pid) {
            _pids.insert(pid);
            _pids.insert(-pid);
            return *this;
        }

    protected:

        /**
         * Method that creates the actual projection for an event.
         */
        void project(const Event& e) {
            // Perform the usual ChargedFinalState projection
            ChargedFinalState::project(e);

            // Loop over the GenParticles and add the ghosts
            std::vector<HepMC::GenParticle *> allParticles = Rivet::particles(e.genEvent());
            foreach(HepMC::GenParticle *g, allParticles) {

                if (_pids.find(g->pdg_id()) != _pids.end()) {
                    Particle p(*g);
                    
                    // A particle with no momentum doesn't have a direction,
                    // set the momentum to a negligibly small value instead.
                    p.setMomentum(1e-6 * p.momentum()/p.momentum().mod());
                    _theParticles.push_back(p);
                }
            }
        }

    private:
        set<PdgId> _pids;
};

} // namespace Rivet

