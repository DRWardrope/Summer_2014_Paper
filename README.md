Summer_2014_Paper
=================
Generation/:    Scripts and configuration files for sample generation.
                Prerequisites: LHAPDF 6.0.6.a0, Fastjet 3.0.6, HepMC 2.06.08, Boost 1.54, Python 2.7.6
                pp->bbbb: Sherpa 2.1.1
                pp->bbcc: Sherpa 2.1.1
                ttbar:    POWHEG-BOX 1 hvq package and Pythia 8.185
Conversion/:    Script to convert Athena pool.root files into HepMC files.
                Requires Athena PhysicsAnalysis/TruthParticleID/McParticleAlgs/ (which uses McParticleTools).
                I used 17.2.11.4
RivetAnalysis/: Analyses HepMC files using Rivet 2.1.2 and produces output ROOT ntuples for use with TMVA and final plotting
RootAnalysis/:  Analyse ROOT files from previous package to produce final results for paper.
                topVeto.cc: Train MVA to reject ttbar events based on output HH signal and ttbar background ROOT files from                                   RivetAnalysis. Reproduces TopVetoWeights\ directory contents.
                FinalPlotter: Applies topVeto MVA and final selection. Produces final plots along with cutflows to std::cout. 
                              Also used to produce TMVA training inputs for final MVA when run with --makeTmvaInput command line                                option.
                TMVAClassification: Acts on output of ./FinalPlotter --makeTmvaInput to train classifier to separate HH from                                          all background. Reproduces weights/ directory contents.
