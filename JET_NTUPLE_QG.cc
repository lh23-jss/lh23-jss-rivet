// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Jet.hh"
#include "TFile.h"
#include "TTree.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <algorithm>

using namespace fastjet;

/// @brief Study of quark and gluon jet substructure in Z+jet and dijet events from pp collisions

namespace Rivet {

  /// @brief Routine for QG substructure analysis
  class JET_NTUPLE_QG : public Analysis {
  public:
    /// Constructor
    JET_NTUPLE_QG() : Analysis("JET_NTUPLE_QG") {}

    /// Book histograms and initialise projections before the run
    void init() {
      _mode = 0;
      if (getOption("MODE") == "DIJET")
        _mode = 1;
      else if (getOption("MODE") == "ZJET")
        _mode = 2;
      else {
        MSG_WARNING("Mode not specified in JET_NTUPLE_QG, using DIJET");
        _mode = 1;
      }

      _jetRadius = getOption("JET_R", 0.4);
      _jetMinPt = getOption("JET_MIN_PT", 50);
      _jetMaxPt = getOption("JET_MAX_PT", 1e9);
      _jetMinAbsEta = getOption("JET_MIN_ABSETA", -1);
      _jetMaxAbsEta = getOption("JET_MAX_ABSETA", 1e9);

      // Initialise and register projections
      FinalState fs(Cuts::abseta < 5 && Cuts::pT > 0 * GeV);
      // Z-jet
      if (_mode == 2) {
        // for the muons
        double mu_pt = 26.;
        double mz_min = (90 - 20);
        double mz_max = (90 + 20);
        double eta_max = 2.4;
        ZFinder zfinder(fs,
                        Cuts::pT > mu_pt * GeV && Cuts::abseta < eta_max,
                        PID::MUON,
                        mz_min * GeV,
                        mz_max * GeV,
                        0.1,
                        ZFinder::ClusterPhotons::NONE,
                        ZFinder::AddPhotons::NO);
        declare(zfinder, "ZFinder");

        eta_max = 2.4;
        FinalState fs_muons(Cuts::abseta < eta_max && Cuts::pT > 0 * GeV);
        IdentifiedFinalState muons_noCut(fs_muons, {PID::MUON, PID::ANTIMUON});
        declare(muons_noCut, "MUONS_NOCUT");

        // Particles for the jets
        VetoedFinalState jet_input(fs);
        jet_input.vetoNeutrinos();
        jet_input.addVetoOnThisFinalState(getProjection<ZFinder>("ZFinder"));
        declare(jet_input, "JET_INPUT");

        FastJets jetfs(jet_input, FastJets::ANTIKT, _jetRadius);
        declare(jetfs, "Jets");
      }
      // dijet
      else {
        // Particles for the jets
        VetoedFinalState jet_input(fs);
        jet_input.vetoNeutrinos();
        declare(jet_input, "JET_INPUT");

        FastJets jetfs(jet_input, FastJets::ANTIKT, _jetRadius);
        declare(jetfs, "Jets");
      }

      // Initialise ROOT file & tree
      _tf = make_unique<TFile>(getOption("ROOTFILE", "Rivet.root").c_str(), "RECREATE");
      _tt = make_unique<TTree>("Rivet", "Rivet physics objects");

      // define branches
      _floatVars["sample"] = 0;
      _floatVars["jet_pt"] = 0;
      _floatVars["jet_eta"] = 0;
      _floatVars["jet_phi"] = 0;
      _floatVars["jet_energy"] = 0;
      _floatVars["jet_nparticles"] = 0;

      for (const auto &v : _lambdaVars) {
        _floatVars[v.name + "_ungroomed"] = 0;
        _floatVars[v.name + "_groomed"] = 0;
      }

      _arrayVars["part_px"];
      _arrayVars["part_py"];
      _arrayVars["part_pz"];
      _arrayVars["part_energy"];
      _arrayVars["part_charge"];
      _arrayVars["part_pid"];

      // book
      for (auto &v : _floatVars) {
        _tt->Branch(v.first.c_str(), &v.second);
      }
      for (auto &v : _arrayVars) {
        _tt->Branch(v.first.c_str(), &v.second, /*bufsize=*/1024000);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event &event) {
      // reset the variable store
      for (auto &v : _floatVars) {
        v.second = 0;
      }
      for (auto &v : _arrayVars) {
        v.second.clear();
      }

      apply<VetoedFinalState>(event, "JET_INPUT");
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 15 * GeV);
      Jets selectedJets;
      JetDefinition jet_def(antikt_algorithm, _jetRadius);

      // z jet
      if (_mode == 2) {
        const FinalState &muons = apply<IdentifiedFinalState>(event, "MUONS_NOCUT");
        if (muons.size() >= 2) {
          Particle muon1 = muons.particlesByPt()[0];
          Particle muon2 = muons.particlesByPt()[1];
          FourMomentum z = muon1.momentum() + muon2.momentum();
        }

        // Reconstruct Z
        const ZFinder &zfinder = apply<ZFinder>(event, "ZFinder");
        if (zfinder.bosons().size() < 1)
          return;

        const Particle &z = zfinder.bosons()[0];
        double zpt = z.pt();

        // Now do selection criteria
        bool passZpJ = false;
        if (jets.size() < 1)
          return;
        const auto &jet1 = jets[0];
        float jet1pt = jet1.pt();
        float asym = fabs((jet1pt - zpt) / (jet1pt + zpt));
        float dphi = Rivet::deltaPhi(jet1.phi(), z.phi());
        passZpJ = ((fabs(jet1.rapidity()) < 1.7) && (zpt > 30) && (asym < 0.3) && (dphi > 2.0));

        if (!passZpJ)
          return;

        // apply jet selection and store it
        if (jet1.pt() < _jetMinPt)
          return;
        if (jet1.pt() > _jetMaxPt)
          return;
        if (jet1.abseta() < _jetMinAbsEta)
          return;
        if (jet1.abseta() > _jetMaxAbsEta)
          return;
        selectedJets.push_back(jet1);
      }
      // di jet
      else {
        bool passDijet = false;
        if (jets.size() < 2)
          return;
        const auto &jet1 = jets.at(0);
        const auto &jet2 = jets.at(1);
        float jet1pt = jet1.pt();
        float jet2pt = jet2.pt();
        float asym = (jet1pt - jet2pt) / (jet1pt + jet2pt);
        float dphi = Rivet::deltaPhi(jet1.phi(), jet2.phi());
        passDijet = ((fabs(jet1.rapidity()) < 1.7) && (fabs(jet2.rapidity()) < 1.7) && (asym < 0.3) && (dphi > 2.0));

        if (!passDijet)
          return;

        // Sort by increasing absolute rapidity
        Jets dijets = {jet1, jet2};
        std::sort(dijets.begin(), dijets.end(), [](const Jet &A, const Jet &B) {
          return fabs(A.rapidity()) < fabs(B.rapidity());
        });

        const auto &centralJet = dijets[0];

        // apply jet selection and store it
        if (centralJet.pt() < _jetMinPt)
          return;
        if (centralJet.pt() > _jetMaxPt)
          return;
        if (centralJet.abseta() < _jetMinAbsEta)
          return;
        if (centralJet.abseta() > _jetMaxAbsEta)
          return;
        selectedJets.push_back(centralJet);
      }

      // save jet constituents
      if (selectedJets.empty())
        return;
      const auto &jet = selectedJets.front();

      _floatVars["sample"] = _mode;
      _floatVars["jet_pt"] = jet.pt();
      _floatVars["jet_eta"] = jet.eta();
      _floatVars["jet_phi"] = jet.phi();
      _floatVars["jet_energy"] = jet.energy();
      _floatVars["jet_nparticles"] = jet.constituents().size();

      // UNGROOMED VERSION
      // -------------------------------------------------------------------
      Particles chargedParticles;
      for (const auto &p : jet.constituents()) {
        if (p.charge() != 0) {
          chargedParticles.push_back(p);
        }
      }
      vector<PseudoJet> chargedJets = jet_def(chargedParticles);

      // Fill hists for each lambda variable
      for (uint lambdaInd = 0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
        const LambdaVar &thisLambdaVar = _lambdaVars[lambdaInd];
        Angularity angularity(thisLambdaVar.beta, _jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
        float val = -1;
        if (thisLambdaVar.isCharged)
          val = (chargedJets.size() > 0) ? angularity(chargedJets[0]) : -1;
        else
          val = angularity(jet);
        _floatVars[thisLambdaVar.name + "_ungroomed"] = val;
      }

      // GROOMED VERSION
      // -------------------------------------------------------------------
      // Get groomed jet
      fastjet::contrib::SoftDrop sd(0, 0.1, _jetRadius);
      PseudoJet groomedJet = sd(jet);
      PseudoJet groomedJetCharged;
      if (chargedJets.size() > 0)
        groomedJetCharged = sd(chargedJets[0]);

      // Fill hists for each lambda variable
      for (uint lambdaInd = 0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
        const LambdaVar &thisLambdaVar = _lambdaVars[lambdaInd];
        Angularity angularity(thisLambdaVar.beta, _jetRadius, thisLambdaVar.kappa, thisLambdaVar.constitCut);
        float val = -1;
        if (thisLambdaVar.isCharged)
          val = (chargedJets.size() > 0) ? angularity(groomedJetCharged) : -1;
        else
          val = angularity(groomedJet);
        _floatVars[thisLambdaVar.name + "_groomed"] = val;
      }

      for (const auto &p : jet.constituents()) {
        _arrayVars["part_px"].push_back(p.px());
        _arrayVars["part_py"].push_back(p.py());
        _arrayVars["part_pz"].push_back(p.pz());
        _arrayVars["part_energy"].push_back(p.energy());
        _arrayVars["part_charge"].push_back(p.charge());
        _arrayVars["part_pid"].push_back(p.pid());
      }

      // fill the tree
      _tt->Fill();

    }  // end analyze() function

    /// Normalise histograms etc., after the run
    void finalize() { _tf->Write(); }  // end of finalize

    /// \class Angularity
    /// definition of angularity
    ///
    class Angularity : public FunctionOfPseudoJet<double> {
    public:
      /// ctor
      Angularity(double alpha, double jet_radius, double kappa = 1.0, Selector constitCut = SelectorPtMin(0.))
          : _alpha(alpha), _radius(jet_radius), _kappa(kappa), _constitCut(constitCut) {}

      /// computation of the angularity itself
      double result(const PseudoJet &jet) const {
        // get the jet constituents
        vector<PseudoJet> constits = jet.constituents();

        // get the reference axis
        PseudoJet reference_axis = _get_reference_axis(jet);

        // do the actual coputation
        double numerator = 0.0, denominator = 0.0;
        for (const auto &c : constits) {
          if (!_constitCut.pass(c))
            continue;
          double pt = c.pt();
          // Note: better compute (dist^2)^(alpha/2) to avoid an extra square root
          numerator += pow(pt, _kappa) * pow(c.squared_distance(reference_axis), 0.5 * _alpha);
          denominator += pt;
        }
        if (denominator == 0)
          return -1;
        // the formula is only correct for the the typical angularities which satisfy either kappa==1 or alpha==0.
        else
          return numerator / (pow(denominator, _kappa) * pow(_radius, _alpha));
      }

    protected:
      PseudoJet _get_reference_axis(const PseudoJet &jet) const {
        if (_alpha > 1)
          return jet;

        Recluster recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme));
        return recluster(jet);
      }

      double _alpha, _radius, _kappa;
      Selector _constitCut;
    };

    /**
         * Lightweight class to hold info about Lambda variable
         */
    class LambdaVar {
    public:
      LambdaVar(const std::string &name_, float kappa_, float beta_, bool isCharged_, Selector constitCut_)
          : name(name_), kappa(kappa_), beta(beta_), isCharged(isCharged_), constitCut(constitCut_) {}

      std::string name;
      float kappa;
      float beta;
      bool isCharged;
      Selector constitCut;
    };

    // This order is important! index in vector used to create YODA plot name
    // Must match that in extracRivetPlotsDijet.py
    const vector<LambdaVar> _lambdaVars = {
        LambdaVar("jet_multiplicity", 0, 0, false, SelectorPtMin(1.)),
        LambdaVar("jet_pTD", 2, 0, false, SelectorPtMin(0.)),
        LambdaVar("jet_LHA", 1, 0.5, false, SelectorPtMin(0.)),
        LambdaVar("jet_width", 1, 1, false, SelectorPtMin(0.)),
        LambdaVar("jet_thrust", 1, 2, false, SelectorPtMin(0.)),
        LambdaVar("jet_multiplicity_charged", 0, 0, true, SelectorPtMin(1.)),
        LambdaVar("jet_pTD_charged", 2, 0, true, SelectorPtMin(0.)),
        LambdaVar("jet_LHA_charged", 1, 0.5, true, SelectorPtMin(0.)),
        LambdaVar("jet_width_charged", 1, 1, true, SelectorPtMin(0.)),
        LambdaVar("jet_thrust_charged", 1, 2, true, SelectorPtMin(0.)),
    };

    // mode for the analysis
    unsigned int _mode;

    // jet radius
    float _jetRadius;
    float _jetMinPt;
    float _jetMaxPt;
    float _jetMinAbsEta;
    float _jetMaxAbsEta;

    /// Output file
    unique_ptr<TFile> _tf;
    /// Output tree
    unique_ptr<TTree> _tt;

    std::map<std::string, float> _floatVars;
    std::map<std::string, std::vector<float>> _arrayVars;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(JET_NTUPLE_QG);

}  // namespace Rivet
