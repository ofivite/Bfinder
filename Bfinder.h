
#ifndef _Bfinder_h
#define _Bfinder_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
// #include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"

double PDG_MUON_MASS       =   0.10565837;
double PDG_PION_MASS       =   0.13957018;
double PDG_PI0_MASS        =   0.1349766;
double PDG_KAON_MASS       =   0.493677;
double PDG_KSHORT_MASS     =   0.497614;
double PDG_KS_MASS         =   0.497614;
double PDG_KSTAR_MASS      =   0.89594;
double PDG_PHI_MASS        =   1.019455;
double PDG_JPSI_MASS       =   3.096916;
double PDG_PSI2S_MASS      =   3.686109;
double PDG_BU_MASS         =   5.27929;
double PDG_B0_MASS         =   5.27961;
double PDG_BS_MASS         =   5.36679;
double PDG_BC_MASS         =   6.2751;
double PDG_LB_MASS         =   5.61951;
double PDG_C               =   29979245800.; // in cm/c
//
// class decleration
//

class Bfinder : public edm::EDAnalyzer {
public:
  explicit Bfinder(const edm::ParameterSet&);
  ~Bfinder();
  int const getMuCat(reco::Muon const& muon) const;

//   void fillPsi(const reco::Candidate& genpsi); // not used
//   void fillV0(const reco::Candidate& genv0); // not used
  bool const HasGoodME11(reco::Muon const& muon, double const dxdzCut) const; // used at some time but not written into tuple


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void CheckL1Triggers(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::string &TrigListNameL1Tmp);
  void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);
  void MatchMuonWithL1L2(const pat::Muon &iMuon, const std::vector<std::string>& TrigListL1L2, std::string &TrigListNameL1L2Tmp);


//   void printout(const RefCountedKinematicVertex& myVertex) const; // not used
//   void printout(const RefCountedKinematicParticle& myParticle) const; // not used
//   void printout(const RefCountedKinematicTree& myTree) const; // not used

//   int PdgIDatTruthLevel(reco::TrackRef Track, edm::Handle<reco::GenParticleCollection> genParticles, int &ParentID); // not used

  float Myctau(const RefCountedKinematicParticle &CandMC, const RefCountedKinematicVertex &DecayVertexMC, //not used
	       const GlobalPoint &PVPtmp, const GlobalError &PVEtmp,float mass_tmp,
	       float &ctau2Dtmp, float &ctauEtemp, float &ctauE2Dtemp );

    // ----------member data ---------------------------
  std::string hlTriggerResults_;
  std::string vtxSample;
  std::string genParticles_;
  std::string muonType;
  std::string muonTypeForPAT;
  bool doMC_;
  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;


  unsigned int             nB;
  unsigned int             nMu;
  unsigned int             nVtx;

  std::vector<std::string>         *triggersMuPL        , *triggersMuML;
  std::vector<std::string>         *triggersL1L2_MuPL   , *triggersL1L2_MuML;

  std::vector<float>    *B0_mass, *pi0_mass;

  std::vector<float>    *B_J_mass;
  std::vector<float>    *B_J_px             , *B_J_py             , *B_J_pz;
  std::vector<float>    *B_J_DecayVtxX      , *B_J_DecayVtxY    , *B_J_DecayVtxZ;
  std::vector<float>    *B_J_DecayVtxXE     , *B_J_DecayVtxYE   , *B_J_DecayVtxZE;

  std::vector<float>    *B_mu_px1_cjp       , *B_mu_py1_cjp     , *B_mu_pz1_cjp;
  std::vector<float>    *B_mu_px2_cjp       , *B_mu_py2_cjp     , *B_mu_pz2_cjp;
  std::vector<float>    *B_J_Prob;

  std::vector<float>    *PV_bestBang_X      , *PV_bestBang_Y    , *PV_bestBang_Z;
  std::vector<float>    *PV_bestBang_XE     , *PV_bestBang_YE   , *PV_bestBang_ZE;
  std::vector<float>    *PV_bestBang_XYE    , *PV_bestBang_XZE  , *PV_bestBang_YZE;
  std::vector<float>    *PV_bestBang_CL;

  std::vector<float>    *PV_bestBang_RF_X   , *PV_bestBang_RF_Y , *PV_bestBang_RF_Z;
  std::vector<float>    *PV_bestBang_RF_XE  , *PV_bestBang_RF_YE, *PV_bestBang_RF_ZE;
  std::vector<float>    *PV_bestBang_RF_XYE , *PV_bestBang_RF_XZE , *PV_bestBang_RF_YZE;
  std::vector<float>    *PV_bestBang_RF_CL;
  std::vector<int>      *PV_bestBang_RF_NTrkDif;

  std::vector<float>    *mum_normChi2       , *mum_dxy_Bsdecay        , *mum_dz_Bsdecay;
  std::vector<int>      *mumCat             , *mumAngT       , *mumNHits  , *mumNPHits;
  std::vector<int>      *mum_isGlobalMuon   , *mum_isTrackerMuon, *mum_isTight,
                        *mum_LastStationOptimizedLowPtT, *mum_LastStationT, *mum_OneStationT,
                        *mum_LastStationAngT, *mum_OneStationAngT, *mum_2DCompatibilityT;
  std::vector<int>      *mum_NMuonHits, *mum_NMuonStations, *mum_NTrackerLayers, *mum_NPixelLayers;
  std::vector<float>    *mum_relIso;

  std::vector<float>    *mup_normChi2       , *mup_dxy_Bsdecay        , *mup_dz_Bsdecay;
  std::vector<int>      *mupCat             , *mupAngT       , *mupNHits  , *mupNPHits;
  std::vector<int>      *mup_isGlobalMuon   , *mup_isTrackerMuon, *mup_isTight,
                        *mup_LastStationOptimizedLowPtT, *mup_LastStationT, *mup_OneStationT,
                        *mup_LastStationAngT, *mup_OneStationAngT, *mup_2DCompatibilityT;
  std::vector<int>      *mup_NMuonHits, *mup_NMuonStations, *mup_NTrackerLayers, *mup_NPixelLayers;
  std::vector<float>    *mup_relIso;

  std::vector<float>    *pho1_Eta, *pho1_Phi, *pho1_Pt, *pho1_Px, *pho1_Py, *pho1_Pz, *pho1_E, *pho1_Et, *pho1_Jpsi_cos;
  std::vector<float>    *pho2_Eta, *pho2_Phi, *pho2_Pt, *pho2_Px, *pho2_Py, *pho2_Pz, *pho2_E, *pho2_Et, *pho2_Jpsi_cos;

///////////////////////

  int                   muAcc, muTrig, weight;

  char triggersL[10000], triggersL1L[10000];

  char triggersMuP[10000],     triggersMuM[10000] ;
  char triggersL1L2_MuP[10000],triggersL1L2_MuM[10000];

  char triggersL1[10000];
  int  nTrgL, nTrgL1L,  nMuonTrgL,  nMuonPTrgL,        nMuonMTrgL;
  int  ntriggersL1L2_MuP, ntriggersL1L2_MuM;


  int  run, event;

};
#endif
