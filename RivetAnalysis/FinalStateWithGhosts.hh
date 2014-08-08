#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

/**
 * FinalState projection that allows injecting "momentumless" particles into
 * the jet clustering. This way you can unambiguously associate hard process
 * quarks to particle jets. Ghost particles are specified by PDG ID.
 */
class FinalStateWithGhosts : public FinalState {

    public:
        
        /**
         * Constructor. Accepts an eta range and a minimal pT.
         */
        FinalStateWithGhosts(double mineta = -100000,
                double maxeta = 100000,
                double minpt = 0.0*GeV) 
            : FinalState(mineta, maxeta, minpt) {}

        /**
         * Clone this projection.
         */
        virtual const Projection *clone() const {
            return new FinalStateWithGhosts(*this);
        }

        /**
         * Add the specified pid and -pid to the list of ghost particles to
         * include.
         */
        FinalStateWithGhosts &ghostIdPair(PdgId pid) {
            _pids.insert(pid);
            _pids.insert(-pid);
            return *this;
        }

    protected:

        /**
         * Method that creates the actual projection for an event.
         */
        void project(const Event& e) {
            // Perform the usual FinalState projection
            FinalState::project(e);

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

