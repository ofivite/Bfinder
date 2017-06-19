// -*- C++ -*-
//
// Package:    Bfinder
// Class:      Bfinder
//
// Author: Jhovanny Andres Mejia,  Eduard De La Cruz Burelo, Ivan Heredia


// system include files
#include <memory>

// user include files
#include "XbFrame/Xb_frame/interface/Bfinder.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"


#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <utility>
#include <string>

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//
// {{{constructors and destructor
//
Bfinder::Bfinder(const edm::ParameterSet& iConfig)
  :
  hlTriggerResults_ (iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT")) ),
  vtxSample ( iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices")) ),
  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  muonType ( iConfig.getUntrackedParameter<std::string>("MuonType",std::string("cleanPatMuons")) ),
  muonTypeForPAT ( iConfig.getUntrackedParameter<std::string>("MuonTypeForPAT",std::string("muons")) ),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0),

  nB(0), nMu(0),  nVtx(0),

   triggersMuPL(0), triggersMuML(0),
   triggersL1L2_MuPL(0), triggersL1L2_MuML(0),

 psiMass(0),

 B_J_mass(0),       B_J_px(0),            B_J_py(0),        B_J_pz(0),
 B_J_DecayVtxX(0),  B_J_DecayVtxY(0),     B_J_DecayVtxZ(0),
 B_J_DecayVtxXE(0), B_J_DecayVtxYE(0),    B_J_DecayVtxZE(0),

 B_mu_px1_cjp(0),   B_mu_py1_cjp(0),  B_mu_pz1_cjp(0),
 B_mu_px2_cjp(0),   B_mu_py2_cjp(0),  B_mu_pz2_cjp(0),

 B_J_Prob(0),

 kaonP_px_0c(0),       kaonP_py_0c(0),  kaonP_pz_0c(0),   kaonP_pt_0c(0),
 kaonP_track_normchi2(0),     kaonP_Hits(0),  kaonP_PHits(0),
 kaonP_NTrackerLayers(0),  kaonP_NPixelLayers(0),

 kaonM_px_0c(0),       kaonM_py_0c(0),  kaonM_pz_0c(0),   kaonM_pt_0c(0),
 kaonM_track_normchi2(0),     kaonM_Hits(0),  kaonM_PHits(0),
 kaonM_NTrackerLayers(0),  kaonM_NPixelLayers(0),

 B_Phi_Prob(0), maxDelta(0),

 PV_bestBang_X(0),      PV_bestBang_Y(0),     PV_bestBang_Z(0),
 PV_bestBang_XE(0),     PV_bestBang_YE(0),    PV_bestBang_ZE(0),
 PV_bestBang_XYE(0),    PV_bestBang_XZE(0),   PV_bestBang_YZE(0),
 PV_bestBang_CL(0),

 PV_bestBang_RF_X(0),   PV_bestBang_RF_Y(0),  PV_bestBang_RF_Z(0),
 PV_bestBang_RF_XE(0),  PV_bestBang_RF_YE(0), PV_bestBang_RF_ZE(0),
 PV_bestBang_RF_XYE(0), PV_bestBang_RF_XZE(0),PV_bestBang_RF_YZE(0),
 PV_bestBang_RF_CL(0),  PV_bestBang_RF_NTrkDif(0),

 mum_normChi2(0),  mum_dxy_Bsdecay(0),    mum_dz_Bsdecay(0)    , mumCat(0) , mumAngT(0)    , mumNHits(0)   , mumNPHits(0),
 mum_isGlobalMuon(0), mum_isTrackerMuon(0), mum_isTight(0),
 mum_LastStationOptimizedLowPtT(0), mum_LastStationT(0), mum_OneStationT(0),
 mum_LastStationAngT(0), mum_OneStationAngT(0), mum_2DCompatibilityT(0),
 mum_NMuonHits(0), mum_NMuonStations(0), mum_NTrackerLayers(0), mum_NPixelLayers(0), mum_relIso(0),

 mup_normChi2(0),  mup_dxy_Bsdecay(0),    mup_dz_Bsdecay(0)    , mupCat(0) , mupAngT(0)    , mupNHits(0)   , mupNPHits(0),
 mup_isGlobalMuon(0), mup_isTrackerMuon(0), mup_isTight(0),
 mup_LastStationOptimizedLowPtT(0), mup_LastStationT(0), mup_OneStationT(0),
 mup_LastStationAngT(0), mup_OneStationAngT(0), mup_2DCompatibilityT(0),
 mup_NMuonHits(0), mup_NMuonStations(0), mup_NTrackerLayers(0), mup_NPixelLayers(0), mup_relIso(0),

 massPiPiPi(0), pfMass(0),
 pionPF_Eta(0),  pionPF_Phi(0),
 pionPF_Pt(0), pi1_maxDeltaR(0), pi2_maxDeltaR(0), pi1_Pt(0),
 pi2_Pt(0), pi1_numberOfGammas(0), pi2_numberOfGammas(0),


 run(0), event(0)

{
   //now do what ever initialization is needed
//   photonCollectionProducer_ = iConfig.getParameter<std::string>("phoProducer");
//   photonCollection_ = iConfig.getParameter<std::string>("photonCollection");
}


Bfinder::~Bfinder()
{

}

//}}}
//
// member functions
// {{{ trigger functions
void Bfinder::CheckHLTTriggers(const std::vector<std::string>& TrigList){ // only do some unuseful shit, define nTRGL and triggersL to write them into tuple

    // using namespace std;
    // using namespace edm;
    // using namespace reco;
    // using namespace helper;


    std::string AllTrg="";
    std::string tmptrig;

    //string alltnames = triggersL;
    //Bool_t checkflag = false;


    int ntrigs=TrigList.size();
    if (ntrigs==0)
        std::cout << "No triggername given in TriggerResults of the input " << std::endl;

    for (int itrig=0; itrig< ntrigs; itrig++) {
        //TString trigName = triggerNames_.triggerName(itrig);
        std::string trigName = TrigList.at(itrig);

	     tmptrig = (std::string) trigName; tmptrig +=" ";
	     AllTrg += tmptrig;
    }

    int m = sprintf(triggersL,"%s","");
    m = sprintf(triggersL,"%s",AllTrg.c_str());
    //cout<<" INFO: Triggers :  "<<m<<"  "<<n<<"  "<<triggersL<<std::endl;
   nTrgL = m;
   //cout<<" INFO: Triggers :  "<<m<<"  "<<nTrgL<<"  "<<triggersL<<std::endl;

   return;
}



void Bfinder::CheckL1Triggers(const edm::Event& iEvent, const edm::EventSetup& iSetup,  std::string &TrigListNameL1Tmp) { // chengest the last argument, then only write the number of trg

  //using namespace std;
   // get L1 trigger info

     edm::ESHandle<L1GtTriggerMenu> menuRcd;
     iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
     const L1GtTriggerMenu* menu = menuRcd.product();

     edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
     iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
     const DecisionWord dWord = gtRecord->decisionWord();

    std::string AllTrgL1="";
    std::string tmptrigL1;

    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {

        std::string trigNameL1 = (algo->second).algoName();

        tmptrigL1 = (std::string) trigNameL1; tmptrigL1 +=" ";
        AllTrgL1 += tmptrigL1;

        //cout << "Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() << std::endl;
    }

    TrigListNameL1Tmp   = AllTrgL1.c_str();

    //  if ( menu->gtAlgorithmResult( "L1_SingleMu3", dWord) )  l1_mu3 = 1;

     return;
}




void Bfinder::MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp){ // changes last agrument to the triggers

  //using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace helper;

  std::string AllTrg="";
  std::string tmptrig;

  int ntrigs=TrigList.size();
  if (ntrigs==0)
    std::cout << "No trigger name given in TriggerResults of the input " << std::endl;

//   cout<<" Trigger in Muon List :";

  for (int itrig=0; itrig< ntrigs; itrig++) {

      std::string trigName = TrigList.at(itrig);

      //if (iMuon.triggerObjectMatchesByPath(trigName.c_str(),(unsigned int)1,(unsigned int)0).empty()==false){
      if (iMuon.triggerObjectMatchesByPath(trigName.c_str(),(unsigned int)0,(unsigned int)0).empty()==false){ //######################################################################################################
	   // cout<<"empty" << endl;
	    //continue;
	// cout<<"In muon: "<<trigName << " " << endl;
      tmptrig = (std::string) trigName; tmptrig +=" ";
      AllTrg += tmptrig;
      }

  }

    //int nMuonTrgLtmp  = AllTrg.size();
    TrigListNameTmp   = AllTrg.c_str();
  //////////////////////////
}


////////////////

void Bfinder::MatchMuonWithL1L2(const pat::Muon &iMuon, const std::vector<std::string>& TrigListL1L2, std::string &TrigListNameL1L2Tmp){


    //using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace helper;


    std::string AllTrgL1L2="";
    std::string tmptrigL1L2;

    int ntrigsL1L2 = TrigListL1L2.size();


    for (int itrigL1L2=0; itrigL1L2< ntrigsL1L2; itrigL1L2++) {

        std::string trigNameL1L2 = TrigListL1L2.at(itrigL1L2);

        //cout<<"i = "<<itrig<<"  "<<trigName.c_str() <<" has " << iMuon.triggerObjectMatchesByPath(trigName.c_str(),(unsigned int)0,(unsigned int)0).empty() << endl;

        if (iMuon.triggerObjectMatchesByFilter(trigNameL1L2.c_str()).empty()==false){
            //continue;
            tmptrigL1L2 = (std::string) trigNameL1L2; tmptrigL1L2 +=" ";
            AllTrgL1L2 += tmptrigL1L2;
        }

    }

    //int nMuonTrgLtmp  = AllTrg.size();
    TrigListNameL1L2Tmp   = AllTrgL1L2.c_str();
    //cout<< "L1/L2: " <<AllTrgL1L2.c_str() << endl;

    //return nMuonTrgLtmp;
    //////////////////////////
}


// }}}



// ------------ method called to for each event  ------------
void Bfinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  //using namespace std;

  //*********************************
  // Get event content information
  //*********************************

  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  Handle< vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);

  Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel(muonType, thePATMuonHandle);

  Handle< reco::PhotonCollection > photonHandle;
  // iEvent.getByLabel("photons", photonHandle);
  iEvent.getByLabel("pfPhotonTranslator", "pfphot", photonHandle);
//  iEvent.getByLabel(photonCollectionProducer_, photonCollection_ , photonHandle);
// const reco::PhotonCollection photonCollection = *(photonHandle.product());

  Handle< vector<reco::PFCandidate> > pfHandle;
  try
{
  iEvent.getByLabel("particleFlowTmp", "AddedMuonsAndHadrons", pfHandle);
}
catch ( ... )
  {
    std::cout << "Couldn't get handle on PFlow!" << std::endl;
  }

  Handle< vector<reco::RecoTauPiZero> > pi0Handle;
  try
{
  iEvent.getByLabel("hpsPFTauProducer", "pizeros", pi0Handle);
}
catch ( ... )
  {
    std::cout << "Couldn't get handle on pi0!" << std::endl;
  }

  //get genParticles
  // Handle<GenParticleCollection> genParticles;
//   if (doMC_)
//     {
  // iEvent.getByLabel(genParticles_, genParticles);
      //    }
    ////////////////////////////////////
 // Get HLT results
  edm::Handle<edm::TriggerResults> hltresults1;
  try {
    std::string const &trig = std::string("TriggerResults::")+hlTriggerResults_;
    iEvent.getByLabel(edm::InputTag(trig),hltresults1);
  }
  catch ( ... )
    {
      std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
    }

    //HLTConfigProvider hltConfig_;

    std::vector<std::string> TrigTable; TrigTable.clear();
    // Get hold of trigger names - based on TriggerResults object
    const edm::TriggerNames& triggerNames1_ = iEvent.triggerNames(*hltresults1);

    for (unsigned int itrig = 0; itrig < hltresults1->size(); ++itrig){
        if ((*hltresults1)[itrig].accept() == 1){
            std::string trigName1 = triggerNames1_.triggerName(itrig);
            //int trigPrescale = hltConfig_.prescaleValue(itrig, trigName1);
            TrigTable.push_back(trigName1);
        }


    }
    CheckHLTTriggers(TrigTable);

    /*
    cout<< "Trigger table: ";
    for( unsigned int i = 0; i<TrigTable.size(); ++i)
        cout<<TrigTable.at(i) << " ";

    cout<< endl;
*/

    // **************************************************
    //trigger condiction "HLT_DoubleMu4_Jpsi_Displaced*"
    // **************************************************

    std::string alltnames = triggersL;
    //Bool_t checkflag = false;
    std::string::size_type trigger1 = alltnames.find("HLT_DoubleMu3p5_LowMass_Displaced",0);
//    std::string::size_type trigger1 = alltnames.find("HLT_DoubleMu4_Jpsi_Displaced",0);
    if(trigger1==std::string::npos)
      {
	//checkflag=true;return;
	return;

      }

    //cout<<"Trigger List "<< triggersL << endl;


    std::vector<std::string> ListTriggerL1L2;

    ListTriggerL1L2.push_back("hltL1MuOpenL1Filtered0");
    ListTriggerL1L2.push_back("hltL2Mu0L2Filtered0");
    ListTriggerL1L2.push_back("hltSingleMu3L2Filtered3");
    ListTriggerL1L2.push_back("hltSingleMu3L3Filtered3");
    ListTriggerL1L2.push_back("hltSingleMu5L3Filtered5");

    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltDiMuonL2PreFiltered0");
    ListTriggerL1L2.push_back("hltDiMuonL3PreFiltered0");
    ListTriggerL1L2.push_back("hltDiMuonL3PreFiltered");
    ListTriggerL1L2.push_back("hltMu0L1MuOpenL3Filtered0");

    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltMu3L1MuOpenL3Filtered3");
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");
    ListTriggerL1L2.push_back("hltMu5L1MuOpenL3Filtered5");
    ListTriggerL1L2.push_back("hltDoubleMuLevel1PathL1OpenFiltered");

    ListTriggerL1L2.push_back("hltL2Mu0L2Filtered0");
    ListTriggerL1L2.push_back("hltMu0TrackJpsiTrackMassFiltered");
    ListTriggerL1L2.push_back("hltMu3TrackJpsiTrackMassFiltered");
    ListTriggerL1L2.push_back("hltMu5TrackJpsiTrackMassFiltered");

    /////////////////////////////////////




  //*********************************
  //Now we get the primary vertex
  //*********************************


  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;
  reco::Vertex BcVtx;

  // get primary vertex
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(vtxSample, recVtxs);
  //unsigned int nVtxTrks = 0;
  bestVtx = *(recVtxs->begin());


  nVtx = recVtxs->size();

  //get primary with beamspot constraint
  Handle<reco::VertexCollection> recVtxes;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", recVtxes);

  bestVtxBS = *(recVtxes->begin());


  //*********************************
  //Now we get the Beam Spot
  //*********************************

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle;
  else std::cout << "No beam spot available from EventSetup" << std::endl;

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  int run1   =  iEvent.id().run();
  int event1 =  iEvent.id().event();


  unsigned int nMu_tmp = thePATMuonHandle->size();

  for( std::vector<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)
    {

      for( std::vector<pat::Muon>::const_iterator iMuon2 = iMuon1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)
	  {
	  if(iMuon1==iMuon2) continue;

	  //opposite charge
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  const pat::Muon *patMuonP = 0;
	  const pat::Muon *patMuonM = 0;
	  TrackRef glbTrackP;
	  TrackRef glbTrackM;
	  TrackRef trkTrackP;
	  TrackRef trkTrackM;

	  if(iMuon1->charge() == 1){ patMuonP = &(*iMuon1); glbTrackP = iMuon1->globalTrack(); trkTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){patMuonM = &(*iMuon1); glbTrackM = iMuon1->globalTrack(); trkTrackM = iMuon1->track();}

	  if(iMuon2->charge() == 1) {patMuonP = &(*iMuon2); glbTrackP = iMuon2->globalTrack(); trkTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){patMuonM = &(*iMuon2); glbTrackM = iMuon2->globalTrack(); trkTrackM = iMuon2->track();}


	  if( trkTrackP.isNull() || trkTrackM.isNull() )
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  TransientTrack muon1TT(trkTrackP, &(*bFieldHandle) );
	  TransientTrack muon2TT(trkTrackM, &(*bFieldHandle) );
      if(!muon1TT.isValid()) continue;
      if(!muon2TT.isValid()) continue;

	  if(iMuon1->pt()<3.5) continue;
	  if(iMuon2->pt()<3.5) continue;

	  if(fabs(iMuon1->eta())>2.2 || fabs(iMuon2->eta())>2.2) continue;

	  if(!(trkTrackM->quality(reco::TrackBase::highPurity))) continue; //quality
	  if(!(trkTrackP->quality(reco::TrackBase::highPurity))) continue; //quality


	  //cout<<"vamos bien 2"<<endl;

      TLorentzVector p4mu1_0c,p4mu2_0c;
      p4mu1_0c.SetPtEtaPhiM(patMuonP->pt(),patMuonP->eta(),patMuonP->phi(), PDG_MUON_MASS);
      p4mu2_0c.SetPtEtaPhiM(patMuonM->pt(),patMuonM->eta(),patMuonM->phi(), PDG_MUON_MASS);
      if (fabs((p4mu1_0c + p4mu2_0c).M() - PDG_JPSI_MASS) > 0.4) continue;


	  //Let's check the vertex and mass
      ParticleMass PM_PDG_MUON_MASS = PDG_MUON_MASS;
//      ParticleMass PM_PDG_JPSI_MASS = PDG_JPSI_MASS;
//      ParticleMass PM_PDG_KAON_MASS = PDG_KAON_MASS;
      ParticleMass PM_PDG_PION_MASS = PDG_PION_MASS;
//      ParticleMass PM_PDG_PI0_MASS  = PDG_PI0_MASS;

      float PM_muon_sigma = PM_PDG_MUON_MASS*1.e-6;
//      float PM_kaon_sigma = PM_PDG_KAON_MASS*1.e-6;
      float PM_pion_sigma = PM_PDG_PION_MASS*1.e-6;
      float chi = 0.;
      float ndf = 0.;

      KinematicParticleFactoryFromTransientTrack pFactory;

      std::vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(pFactory.particle(muon1TT, PM_PDG_MUON_MASS, chi,ndf, PM_muon_sigma));
      muonParticles.push_back(pFactory.particle(muon2TT, PM_PDG_MUON_MASS, chi,ndf, PM_muon_sigma));

      KinematicParticleVertexFitter   fitter;
      RefCountedKinematicTree         ujVFT;
      ujVFT = fitter.fit(muonParticles);
      if (!ujVFT->isValid()) continue;

      ujVFT->movePointerToTheTop();
      RefCountedKinematicParticle MUMUparticle       = ujVFT->currentParticle();
      RefCountedKinematicVertex   MUMUvtx            = ujVFT->currentDecayVertex();
      ujVFT->movePointerToTheFirstChild();
      RefCountedKinematicParticle MUP_cMUMU          = ujVFT->currentParticle();
      ujVFT->movePointerToTheNextChild();
      RefCountedKinematicParticle MUM_cMUMU          = ujVFT->currentParticle();
      if (!MUMUvtx->vertexIsValid())  continue;
      double JP_Prob_tmp = TMath::Prob(MUMUvtx->chiSquared(), MUMUvtx->degreesOfFreedom());
      if(JP_Prob_tmp < 0.01) continue;

      double MUMU_mass_c0 = MUMUparticle->currentState().mass();
      if ( MUMU_mass_c0 < PDG_JPSI_MASS - 0.15 ) continue;
      if ( MUMU_mass_c0 > PDG_JPSI_MASS + 0.15 ) continue;


      // Get Matching Muon to HLT
      std::string ListTriggMuP_tmp="";
      std::string ListTriggMuM_tmp="";

      if(iMuon1->charge()== 1) MatchMuonWithTriggers(*iMuon1, TrigTable, ListTriggMuP_tmp);
      if(iMuon1->charge()==-1) MatchMuonWithTriggers(*iMuon1, TrigTable, ListTriggMuM_tmp);
      if(iMuon2->charge()== 1) MatchMuonWithTriggers(*iMuon2, TrigTable, ListTriggMuP_tmp);
      if(iMuon2->charge()==-1) MatchMuonWithTriggers(*iMuon2, TrigTable, ListTriggMuM_tmp);

      // *** *************** MUON MACHINT HERE ********** ****************************

      std::string::size_type triggerMupOS = ListTriggMuP_tmp.find("HLT_DoubleMu3p5_LowMass_Displaced",0);
      std::string::size_type triggerMuNeg = ListTriggMuM_tmp.find("HLT_DoubleMu3p5_LowMass_Displaced",0);

      if(triggerMupOS==std::string::npos || triggerMuNeg==std::string::npos )
        {
          continue;
        }

      // ******************************* end muon machint *******************************

      int nMuonP = sprintf(triggersMuP,"%s",ListTriggMuP_tmp.c_str());
      int nMuonM = sprintf(triggersMuM,"%s",ListTriggMuM_tmp.c_str());

      nMuonPTrgL = nMuonP;
      nMuonMTrgL = nMuonM;


/*
RefCountedKinematicTree

    ->currentParticle();
    RefCountedKinematicParticle
    https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_6/doc/html/d2/dff/classKinematicParticle.html

        -> currentState()
        KinematicState
        https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_6/doc/html/db/d18/classKinematicState.html
            .globalMomentum().x()

            -> mass()
            -> kinematicParameters ()
            https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_6/doc/html/d6/dc6/classKinematicParameters.html
                .mass()

                .momentum()
                https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_6/doc/html/de/d6a/classVector3DBase.html
                    . x()
                    . y()

            -> kinematicParametersError ()
            https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_6/doc/html/d4/de9/classKinematicParametersError.html#a37281e3886c52d6fa55eea05218cd583

    ->currentDecayVertex()
    RefCountedKinematicVertex
    https://cmssdt.cern.ch/SDT/doxygen/CMSSW_2_2_13_HLT/doc/html/df/dd9/classKinematicTree.html

        ->vertexIsValid()

        ->position().
            ->x()
*/
// 	   int PId1=0; int PId2=0;




	   pat::GenericParticle patTrack_Kp;
 	   pat::GenericParticle patTrack_Km;

	   for(vector<pat::GenericParticle>::const_iterator iKaonP = thePATTrackHandle->begin(); iKaonP != thePATTrackHandle->end(); ++iKaonP )
	     {

// 	       int ngenT1 = 0;//PdgIDatTruthLevel(iKaonP->track(), genParticles, PId1);
	       patTrack_Kp = *iKaonP;

        if(iKaonP->pt() < 2. || iKaonP->charge() != 1 || fabs(iKaonP->eta()) > 2.5) continue;


	       if(!(patTrack_Kp.track()->quality(reco::TrackBase::highPurity))) continue;

               bool matchflag = false;
               const reco::CandidatePtrVector & mu1P_overlaps = patTrack_Kp.overlaps(muonTypeForPAT);
               if ( mu1P_overlaps.size() > 0 ) //std::cout << "patTrack_Kp overlaps with a muon." << endl;
               for (size_t i = 0; i < mu1P_overlaps.size(); ++i) {
                 const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu1P_overlaps[i]);
                 if (mu) {
                   // check here that muon match isn't the same as a muon used in the reco...
                   if (mu==patMuonP || mu==patMuonM)  {
                       //std::cout << "match between patTrack_Kp and patMuonP/patMuonM" << endl;
                       matchflag=true;
                   }
                 }
               }

               if(matchflag) continue;

               TransientTrack kaonPTT(patTrack_Kp.track(), &(*bFieldHandle) );
               if(!kaonPTT.isValid()) continue;

/////
for(vector<pat::GenericParticle>::const_iterator iKaonM = thePATTrackHandle->begin(); iKaonM != thePATTrackHandle->end(); ++iKaonM )
  {
   if(iKaonP == iKaonM) continue;
// 	       int ngenT1 = 0;//PdgIDatTruthLevel(iKaonP->track(), genParticles, PId1);
   patTrack_Km = *iKaonM;
   if(iKaonM->pt() < 2. || iKaonM->charge() != -1 || fabs(iKaonP->eta()) > 2.5) continue;


    if(!(patTrack_Km.track()->quality(reco::TrackBase::highPurity))) continue;

          bool matchflag = false;
          const reco::CandidatePtrVector & mu2P_overlaps = patTrack_Km.overlaps(muonTypeForPAT);
          if ( mu2P_overlaps.size() > 0 ) //std::cout << "patTrack_Kp overlaps with a muon." << endl;
          for (size_t i = 0; i < mu2P_overlaps.size(); ++i) {
            const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu2P_overlaps[i]);
            if (mu) {
              // check here that muon match isn't the same as a muon used in the reco...
              if (mu==patMuonP || mu==patMuonM)  {
                  //std::cout << "match between patTrack_Kp and patMuonP/patMuonM" << endl;
                  matchflag=true;
              }
            }
          }

          if(matchflag) continue;

          TransientTrack kaonMTT(patTrack_Km.track(), &(*bFieldHandle) );
          if(!kaonMTT.isValid()) continue;


               TLorentzVector p4kaonP_0c, p4kaonM_0c, p4phi_0c, p4jpsi_0c;
               p4kaonP_0c.SetXYZM(iKaonP->px(),iKaonP->py(),iKaonP->pz(), PDG_PION_MASS);
               p4kaonM_0c.SetXYZM(iKaonM->px(),iKaonM->py(),iKaonM->pz(), PDG_PION_MASS);
	       p4jpsi_0c = p4mu1_0c + p4mu2_0c;

                std::vector<RefCountedKinematicParticle> phiParticles;
                phiParticles.push_back(pFactory.particle(kaonPTT, PM_PDG_PION_MASS, chi,ndf, PM_pion_sigma));
                phiParticles.push_back(pFactory.particle(kaonMTT, PM_PDG_PION_MASS, chi,ndf, PM_pion_sigma));

                KinematicParticleVertexFitter   phifitter;
                RefCountedKinematicTree         phiTree;
                phiTree = phifitter.fit(phiParticles);
                if (!phiTree->isValid()) continue;

                phiTree->movePointerToTheTop();
                RefCountedKinematicParticle PHIparticle       = phiTree->currentParticle();
                RefCountedKinematicVertex   PHIvtx            = phiTree->currentDecayVertex();
                if (!PHIvtx->vertexIsValid())  continue;
                double phi_Prob_tmp = TMath::Prob(PHIvtx->chiSquared(), PHIvtx->degreesOfFreedom());
                if(phi_Prob_tmp < 0.1) continue;

//                double PHI_mass_c0 = PHIparticle->currentState().mass();



         for( std::vector<reco::RecoTauPiZero>::const_iterator iPiZero1 = pi0Handle->begin(); iPiZero1 != pi0Handle->end(); ++iPiZero1)
         {
  //            if (iPiZero1->particleId() != 5) continue;

              // reco::RecoTauPiZero localpionPF = reco::PFCandidate(*iPiZero1);
              TLorentzVector PiZeroP4_1, PiZeroP4_2, psiP4;
              if ( iPiZero1->pt() < 2. || iPiZero1->numberOfGammas() != 2) continue;

              PiZeroP4_1.SetPxPyPzE(iPiZero1->px(), iPiZero1->py(), iPiZero1->pz(), iPiZero1->energy());
              if ( (PiZeroP4_1 + p4kaonP_0c + p4kaonM_0c).M() < 0.3 || (PiZeroP4_1 + p4kaonP_0c + p4kaonM_0c).M() > 1.) continue;

//              PiZeroP4_2.SetPxPyPzE(iPiZero2->px(), iPiZero2->py(), iPiZero2->pz(), iPiZero2->energy());
//              psiP4 = PiZeroP4_1 + PiZeroP4_2 + p4mu1_0c + p4mu2_0c;
//              if ( fabs(psiP4.M() - PDG_PSI2S_MASS) > 0.1 ) continue;


              pfMass->push_back(iPiZero1->mass());

                // pionPF_Eta->push_back( localpionPF.eta() );
                // pionPF_Phi->push_back( localpionPF.phi() );
                // pionPF_Pt->push_back( localpionPF.pt() );

              massPiPiPi->push_back( (PiZeroP4_1 + p4kaonP_0c + p4kaonM_0c).M()  );
//              psiMass->push_back( psiP4.M() );
              pi1_maxDeltaR->push_back( sqrt(iPiZero1->maxDeltaEta()*iPiZero1->maxDeltaEta() + iPiZero1->maxDeltaPhi()*iPiZero1->maxDeltaPhi()) );
//              pi2_maxDeltaR->push_back( sqrt(iPiZero2->maxDeltaEta()*iPiZero2->maxDeltaEta() + iPiZero2->maxDeltaPhi()*iPiZero2->maxDeltaPhi()) );
              pi1_Pt->push_back( PiZeroP4_1.Pt() );
//              pi2_Pt->push_back( localpionPF.energy() );
              pi1_numberOfGammas->push_back( iPiZero1->numberOfGammas() );

//              pi2_numberOfGammas->push_back( iPiZero2->numberOfGammas() );


              reco::Vertex bestPV_Bang;

              Double_t PV_bestBang_X_temp  = -10000.0;
              Double_t PV_bestBang_Y_temp  = -10000.0;
              Double_t PV_bestBang_Z_temp  = -10000.0;
              Double_t PV_bestBang_XE_temp = -10000.0;
              Double_t PV_bestBang_YE_temp = -10000.0;
              Double_t PV_bestBang_ZE_temp = -10000.0;
              Double_t PV_bestBang_XYE_temp= -10000.0;
              Double_t PV_bestBang_XZE_temp= -10000.0;
              Double_t PV_bestBang_YZE_temp= -10000.0;
              Double_t PV_bestBang_CL_temp = -10000.0;

              Double_t lip = -100000.0;

              for(size_t i = 0; i < recVtxes->size(); ++i)
              {
                   const Vertex &PVtxBeSp = (*recVtxes)[i];

                   Double_t dx = (*MUMUvtx).position().x() - PVtxBeSp.x();
                   Double_t dy = (*MUMUvtx).position().y() - PVtxBeSp.y();
                   Double_t dz = (*MUMUvtx).position().z() - PVtxBeSp.z();
                   Double_t cosAlphaXYZ = ( p4jpsi_0c.Px() * dx + p4jpsi_0c.Py()*dy + p4jpsi_0c.Pz()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* p4jpsi_0c.P() );

                   if(cosAlphaXYZ>lip)
                   {
                       lip = cosAlphaXYZ ;

                       PV_bestBang_X_temp     = PVtxBeSp.x();
                       PV_bestBang_Y_temp     = PVtxBeSp.y();
                       PV_bestBang_Z_temp     = PVtxBeSp.z();
                       PV_bestBang_XE_temp    = PVtxBeSp.covariance(0, 0);
                       PV_bestBang_YE_temp    = PVtxBeSp.covariance(1, 1);
                       PV_bestBang_ZE_temp    = PVtxBeSp.covariance(2, 2);
                       PV_bestBang_XYE_temp   = PVtxBeSp.covariance(0, 1);
                       PV_bestBang_XZE_temp   = PVtxBeSp.covariance(0, 2);
                       PV_bestBang_YZE_temp   = PVtxBeSp.covariance(1, 2);
                       PV_bestBang_CL_temp    = (TMath::Prob(PVtxBeSp.chi2(),(int)PVtxBeSp.ndof()) );

                       bestPV_Bang = PVtxBeSp;
                   }
              }

	    reco::Vertex bestVtxRf = bestPV_Bang;
/*
             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
             /////////////////////////////////// // try refitting the primary without the tracks in the B reco candidate

              // first get tracks from the original primary
              vector<reco::TransientTrack> vertexTracks;

               for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestPV_Bang.tracks_begin();
                 iTrack != bestPV_Bang.tracks_end(); ++iTrack)
              {
                    // compare primary tracks to check for matches with B cand
                    TrackRef trackRef = iTrack->castTo<TrackRef>();

                    // the 4 tracks in the B cand are  patTrack_Kp glbTrackP glbTrackM
                    if (  !(
                               (trkTrackP==trackRef)            ||
                               (trkTrackM==trackRef)           ) )
                       {
                           TransientTrack tt(trackRef, &(*bFieldHandle) );
                           vertexTracks.push_back(tt);
                       } //else { std::cout << "found track match with primary" << endl;}
               }

               // if no tracks in primary or no reco track included in primary then don't do anything
               // if so, then update bctau_temp and bctauMPV_temp

               reco::Vertex bestVtxRf = bestPV_Bang;
               GlobalPoint PVRfP = GlobalPoint( bestPV_Bang.x(), bestPV_Bang.y(), bestPV_Bang.z() );

               if (  vertexTracks.size()>0 && (bestPV_Bang.tracksSize()!=vertexTracks.size()) ) {
                 AdaptiveVertexFitter theFitter;

                 TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);

                    if ( v.isValid() ) {

                  //calculate ctau with the new vertex to compare to the old one.
                  //GlobalPoint PVRfP = GlobalPoint( v.position().x(), v.position().y(), v.position().z() );
                  //reco::Vertex recoV = (reco::Vertex)v;

                  //GlobalError PVRfE = GlobalError( recoV.error() );
                  //bctauRf_temp = Myctau(bCandCjp, BsDecayVertexCjp, PVRfP, PVRfE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);

                  //set bestVtxRf as new best vertex to fill variables for ntuple
                  bestVtxRf = reco::Vertex(v);
                }
               }

*/

               //cout << "PV bestPV_Bang: " <<bestPV_Bang.x()<< " "<<bestPV_Bang.y()<<" "<<bestPV_Bang.z()<< endl;



               if(iMuon1->charge() == 1 ) mupCategory = getMuCat( *iMuon1 );
               if(iMuon1->charge() == -1) mumCategory = getMuCat( *iMuon1 );
               if(iMuon2->charge() == 1 ) mupCategory = getMuCat( *iMuon2 );
               if(iMuon2->charge() == -1) mumCategory = getMuCat( *iMuon2 );

               const reco::Muon *recoMuonM = patMuonM;
               const reco::Muon *recoMuonP = patMuonP;

               // fill candidate variables now
               //            PVindex_temp = PV_bestBang_X->size();
                          PV_bestBang_X    ->push_back(  PV_bestBang_X_temp);
                          PV_bestBang_Y    ->push_back(  PV_bestBang_Y_temp);
                          PV_bestBang_Z    ->push_back(  PV_bestBang_Z_temp);
                          PV_bestBang_XE   ->push_back(  PV_bestBang_XE_temp);
                          PV_bestBang_YE   ->push_back(  PV_bestBang_YE_temp);
                          PV_bestBang_ZE   ->push_back(  PV_bestBang_ZE_temp);
                          PV_bestBang_XYE  ->push_back(  PV_bestBang_XYE_temp);
                          PV_bestBang_XZE  ->push_back(  PV_bestBang_XZE_temp);
                          PV_bestBang_YZE  ->push_back(  PV_bestBang_YZE_temp);
                          PV_bestBang_CL   ->push_back(  PV_bestBang_CL_temp);
                          //
                          PV_bestBang_RF_X ->push_back(    bestVtxRf.x() );
                          PV_bestBang_RF_Y ->push_back(    bestVtxRf.y() );
                          PV_bestBang_RF_Z ->push_back(    bestVtxRf.z() );
                          PV_bestBang_RF_XE->push_back(    bestVtxRf.covariance(0, 0) );
                          PV_bestBang_RF_YE->push_back(    bestVtxRf.covariance(1, 1) );
                          PV_bestBang_RF_ZE->push_back(    bestVtxRf.covariance(2, 2) );
                          PV_bestBang_RF_XYE->push_back(   bestVtxRf.covariance(0, 1) );
                          PV_bestBang_RF_XZE->push_back(   bestVtxRf.covariance(0, 2) );
                          PV_bestBang_RF_YZE->push_back(   bestVtxRf.covariance(1, 2) );
                          PV_bestBang_RF_CL->push_back(    ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );

               // Only save the first time
               if(nB==0){

                 // Get L1 trigger to level Event
                 std::string ListTriggL1_tmp="";
                 CheckL1Triggers(iEvent, iSetup, ListTriggL1_tmp);

                 nTrgL1L = ListTriggL1_tmp.size();
                 nMu  = nMu_tmp;
                 // cout<< "*Number of Muons : " << nMu_tmp << endl;

                 run   =  run1;
                 event =  event1;

               } // end nB==0


               // Get List for L1/L2 triggers matching to Muon
               std::string ListTriggL1L2_MuP="";
               std::string ListTriggL1L2_MuM="";
               if(iMuon1->charge()== 1)  MatchMuonWithL1L2(*iMuon1, ListTriggerL1L2, ListTriggL1L2_MuP);
               if(iMuon1->charge()== -1) MatchMuonWithL1L2(*iMuon1, ListTriggerL1L2, ListTriggL1L2_MuM);
               if(iMuon2->charge()== 1)  MatchMuonWithL1L2(*iMuon2, ListTriggerL1L2, ListTriggL1L2_MuP);
               if(iMuon2->charge()== -1) MatchMuonWithL1L2(*iMuon2, ListTriggerL1L2, ListTriggL1L2_MuM);

               int nL1L2MuP = sprintf(triggersL1L2_MuP,"%s",ListTriggL1L2_MuP.c_str());
               int nL1L2MuM = sprintf(triggersL1L2_MuM,"%s",ListTriggL1L2_MuM.c_str());

               triggersL1L2_MuPL->push_back(triggersL1L2_MuP);
               triggersL1L2_MuML->push_back(triggersL1L2_MuM);

               ntriggersL1L2_MuP = nL1L2MuP;
               ntriggersL1L2_MuM = nL1L2MuM;


               triggersMuPL ->push_back(triggersMuP);
               triggersMuML ->push_back(triggersMuM);


 //              PV_bestBang_RF_NTrkDif->push_back( bestPV_Bang.tracksSize() - vertexTracks.size() );
               /////////////////////////////////////////////////////////////////////////////////////////////
               ////////////////////////////////////////////// HERE INDEX IS [nB] /////////////////////////
                        //


                          B_J_mass         ->push_back( MUMU_mass_c0 );
                          B_J_px           ->push_back( MUMUparticle->currentState().globalMomentum().x() );
                          B_J_py           ->push_back( MUMUparticle->currentState().globalMomentum().y() );
                          B_J_pz           ->push_back( MUMUparticle->currentState().globalMomentum().z() );
                          B_J_DecayVtxX    ->push_back( MUMUvtx->position().x() );
                          B_J_DecayVtxY    ->push_back( MUMUvtx->position().y() );
                          B_J_DecayVtxZ    ->push_back( MUMUvtx->position().z() );
                          B_J_DecayVtxXE   ->push_back( MUMUvtx->error().cxx() );
                          B_J_DecayVtxYE   ->push_back( MUMUvtx->error().cyy() );
                          B_J_DecayVtxZE   ->push_back( MUMUvtx->error().czz() );

                          B_mu_px1_cjp ->push_back(MUP_cMUMU->currentState().kinematicParameters().momentum().x());
                          B_mu_py1_cjp ->push_back(MUP_cMUMU->currentState().kinematicParameters().momentum().y());
                          B_mu_pz1_cjp ->push_back(MUP_cMUMU->currentState().kinematicParameters().momentum().z());
           //                B_J_charge1->push_back(mНУ А как же всем известный АМХ 40 ?)))u1CandMC->currentState().particleCharge()); // noneed

                          B_mu_px2_cjp ->push_back(MUM_cMUMU->currentState().kinematicParameters().momentum().x());
                          B_mu_py2_cjp ->push_back(MUM_cMUMU->currentState().kinematicParameters().momentum().y());
                          B_mu_pz2_cjp ->push_back(MUM_cMUMU->currentState().kinematicParameters().momentum().z());
           //                B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

                          B_J_Prob  ->push_back(JP_Prob_tmp);


   //------------------//

                          kaonP_px_0c       ->push_back(iKaonP->px());
                          kaonP_py_0c       ->push_back(iKaonP->py());
                          kaonP_pz_0c       ->push_back(iKaonP->pz());
                          kaonP_pt_0c       ->push_back(iKaonP->pt());
	    
                          kaonP_track_normchi2  ->push_back(patTrack_Kp.track()->normalizedChi2());
                          kaonP_Hits       ->push_back(patTrack_Kp.track()->numberOfValidHits() );
                          kaonP_PHits      ->push_back(patTrack_Kp.track()->hitPattern().numberOfValidPixelHits() );
	                  kaonP_NTrackerLayers->push_back ( patTrack_Kp.track()->hitPattern().trackerLayersWithMeasurement() );
   		          kaonP_NPixelLayers->push_back ( patTrack_Kp.track()->hitPattern().pixelLayersWithMeasurement() );

                          kaonM_px_0c       ->push_back(iKaonM->px());
                          kaonM_py_0c       ->push_back(iKaonM->py());
                          kaonM_pz_0c       ->push_back(iKaonM->pz());
                          kaonM_pt_0c       ->push_back(iKaonM->pt());

                          kaonM_track_normchi2  ->push_back(patTrack_Km.track()->normalizedChi2());
                          kaonM_Hits       ->push_back(patTrack_Km.track()->numberOfValidHits() );
                          kaonM_PHits      ->push_back(patTrack_Km.track()->hitPattern().numberOfValidPixelHits() );
	 		  kaonM_NTrackerLayers->push_back ( patTrack_Km.track()->hitPattern().trackerLayersWithMeasurement() );
   		          kaonM_NPixelLayers->push_back ( patTrack_Km.track()->hitPattern().pixelLayersWithMeasurement() );
		
			  B_Phi_Prob->push_back(phi_Prob_tmp);
			  maxDelta->push_back(std::max( std::max(PiZeroP4_1.DeltaR(p4kaonP_0c), PiZeroP4_1.DeltaR(p4kaonM_0c)), p4kaonP_0c.DeltaR(p4kaonM_0c) ));
   //------------------//
                  mumCat->push_back( mumCategory );
   	              mum_isGlobalMuon->push_back ( recoMuonM->isGlobalMuon() );
                  mum_isTrackerMuon->push_back ( recoMuonM->isTrackerMuon() );
                  mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
/*                  if ( BsTV.isValid() )
        	          mum_isTight->push_back ( muon::isTightMuon(*recoMuonM, BsVtx ) );
                  else
*/                    mum_isTight->push_back ( -1 );
//                  mum_dxy_Bsdecay->push_back( recoMuonM->muonBestTrack()->dxy(BsDecayPoint) );// el dxy del Muon negatico respetcto del PV con BSc (el de mayor pt)
//                  mum_dz_Bsdecay->push_back( recoMuonM->muonBestTrack()->dz(BsDecayPoint) );

         	       if (!(glbTrackM.isNull()))
                    {
   	                mum_normChi2->push_back( glbTrackM->normalizedChi2() );
   		              mum_NMuonHits->push_back ( glbTrackM->hitPattern().numberOfValidMuonHits() );
   		             }
   	             else
   		             {
   	 	                mum_normChi2->push_back(-1);
   		                mum_NMuonHits->push_back(-1);
   		              }

               		mum_NMuonStations->push_back ( recoMuonM->numberOfMatchedStations() );

               		mum_NTrackerLayers->push_back ( trkTrackM->hitPattern().trackerLayersWithMeasurement() );
               		mum_NPixelLayers->push_back ( trkTrackM->hitPattern().pixelLayersWithMeasurement() );
                   mumNPHits->push_back( trkTrackM->hitPattern().numberOfValidPixelHits() );
                   mumNHits->push_back( trkTrackM->numberOfValidHits() );

               		float mum_chIso  = patMuonM->chargedHadronIso();
               		float mum_nhIso  = patMuonM->neutralHadronIso();
               		float mum_phIso  = patMuonM->photonIso();
               		float mum_puIso  = patMuonM->puChargedHadronIso();
               		mum_relIso->push_back( (mum_chIso + std::max(0.0, mum_nhIso + mum_phIso - 0.5 * mum_puIso)) / patMuonM->pt() );

                  mum_LastStationOptimizedLowPtT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMLastStationOptimizedLowPtTight) );
                  mum_LastStationT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMLastStationTight) );
                  mum_OneStationT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationTight) );
                  mum_LastStationAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMLastStationAngTight) );
                  mum_OneStationAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationAngTight) );
                  mum_2DCompatibilityT->push_back( muon::isGoodMuon(*recoMuonM,muon::TM2DCompatibilityTight) );




   //------------------//
                  mupCat->push_back( mupCategory );
                  mup_isGlobalMuon->push_back ( recoMuonP->isGlobalMuon() );
                  mup_isTrackerMuon->push_back ( recoMuonP->isTrackerMuon() );
                  mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
/*                  if ( BsTV.isValid() )
        	           mup_isTight->push_back ( muon::isTightMuon(*recoMuonP, BsVtx ) );
                  else
*/                    mup_isTight->push_back ( -1 );

//                  mup_dxy_Bsdecay->push_back( recoMuonP->muonBestTrack()->dxy(BsDecayPoint) );// el dxy del Muon negatico respetcto del PV con BSc (el de mayor pt)
//                  mup_dz_Bsdecay->push_back( recoMuonP->muonBestTrack()->dz(BsDecayPoint) );

         	       if (!(glbTrackP.isNull()))
                  {
   	                mup_normChi2->push_back( glbTrackP->normalizedChi2() );
   		              mup_NMuonHits->push_back ( glbTrackP->hitPattern().numberOfValidMuonHits() );
   		            }
   	             else
   		            {
   	 	              mup_normChi2->push_back(-1);
   		              mup_NMuonHits->push_back (-1);
   		            }
   		              mup_NMuonStations->push_back ( recoMuonP->numberOfMatchedStations() );

               		mup_NTrackerLayers->push_back ( trkTrackP->hitPattern().trackerLayersWithMeasurement() );
               		mup_NPixelLayers->push_back ( trkTrackP->hitPattern().pixelLayersWithMeasurement() );
                   mupNPHits->push_back( trkTrackP->hitPattern().numberOfValidPixelHits() );
                   mupNHits->push_back( trkTrackP->numberOfValidHits() );

               		float mup_chIso  = patMuonP->chargedHadronIso();
               		float mup_nhIso  = patMuonP->neutralHadronIso();
               		float mup_phIso  = patMuonP->photonIso();
               		float mup_puIso  = patMuonP->puChargedHadronIso();
               		mup_relIso->push_back( (mup_chIso + std::max(0.0, mup_nhIso + mup_phIso - 0.5 * mup_puIso)) / patMuonP->pt() );

                  mup_LastStationOptimizedLowPtT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMLastStationOptimizedLowPtTight) );
                  mup_LastStationT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMLastStationTight) );
                  mup_OneStationT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationTight) );
                  mup_LastStationAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMLastStationAngTight) );
                  mup_OneStationAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationAngTight) );
                  mup_2DCompatibilityT->push_back( muon::isGoodMuon(*recoMuonP,muon::TM2DCompatibilityTight) );


                          nB++;

                          /////////////////////////////////////////////////

                          muonParticles.clear();
 //                         vertexTracks.clear();
                  } // pion
              } // pion
	} // pion
   } // muon from jpsi
   }//muon from jpsi


      //fill the tree and clear the vectors
      if (nB > 0 )
        {

          //std::cout << "filling tree" << endl;
          tree_->Fill();
        }

      // {{{ // clearence of the variables

      nB = 0;  nMu = 0;    nVtx = 0;

      //triggersL1L->clear();
      triggersMuPL->clear(); triggersMuML->clear();
      triggersL1L2_MuPL->clear(); triggersL1L2_MuML->clear();

      psiMass->clear();

      B_J_mass->clear();       B_J_px->clear();            B_J_py->clear();        B_J_pz->clear();
      B_J_DecayVtxX->clear();  B_J_DecayVtxY->clear();     B_J_DecayVtxZ->clear();
      B_J_DecayVtxXE->clear(); B_J_DecayVtxYE->clear();    B_J_DecayVtxZE->clear();


      B_mu_px1_cjp->clear();   B_mu_py1_cjp->clear();  B_mu_pz1_cjp->clear();
      B_mu_px2_cjp->clear();   B_mu_py2_cjp->clear();  B_mu_pz2_cjp->clear();

      B_J_Prob->clear();
      //

      kaonP_px_0c->clear();     kaonP_py_0c->clear();    kaonP_pz_0c->clear();   kaonP_pt_0c->clear();
      kaonP_track_normchi2->clear();   kaonP_Hits->clear();    kaonP_PHits->clear();
      kaonP_NTrackerLayers->clear();  kaonP_NPixelLayers->clear();

      kaonM_px_0c->clear();     kaonM_py_0c->clear();    kaonM_pz_0c->clear();   kaonM_pt_0c->clear();
      kaonM_track_normchi2->clear();   kaonM_Hits->clear();    kaonM_PHits->clear();
      kaonM_NTrackerLayers->clear();  kaonM_NPixelLayers->clear();

      B_Phi_Prob->clear(); maxDelta->clear();

      PV_bestBang_X->clear();      PV_bestBang_Y->clear();     PV_bestBang_Z->clear();
      PV_bestBang_XE->clear();     PV_bestBang_YE->clear();    PV_bestBang_ZE->clear();
      PV_bestBang_XYE->clear();    PV_bestBang_XZE->clear();   PV_bestBang_YZE->clear();
      PV_bestBang_CL->clear();

      PV_bestBang_RF_X->clear();   PV_bestBang_RF_Y->clear();  PV_bestBang_RF_Z->clear();
      PV_bestBang_RF_XE->clear();  PV_bestBang_RF_YE->clear(); PV_bestBang_RF_ZE->clear();
      PV_bestBang_RF_XYE->clear(); PV_bestBang_RF_XZE->clear();PV_bestBang_RF_YZE->clear();
      PV_bestBang_RF_CL->clear();  PV_bestBang_RF_NTrkDif->clear();

      mum_normChi2->clear();  mum_dxy_Bsdecay->clear();    mum_dz_Bsdecay->clear(); mumCat->clear();    mumAngT->clear();   mumNHits->clear();  mumNPHits->clear();
      mum_isGlobalMuon->clear(); mum_isTrackerMuon->clear(); mum_isTight->clear();
      mum_LastStationOptimizedLowPtT->clear(); mum_LastStationT->clear(); mum_OneStationT->clear();
      mum_LastStationAngT->clear(); mum_OneStationAngT->clear(); mum_2DCompatibilityT->clear();
      mum_NMuonHits->clear(); mum_NMuonStations->clear(); mum_NTrackerLayers->clear(); mum_NPixelLayers->clear(); mum_relIso->clear();

      mup_normChi2->clear();   mup_dxy_Bsdecay->clear();    mup_dz_Bsdecay->clear(); mupCat->clear();    mupAngT->clear();   mupNHits->clear();  mupNPHits->clear();
      mup_isGlobalMuon->clear(); mup_isTrackerMuon->clear(); mup_isTight->clear();
      mup_LastStationOptimizedLowPtT->clear(); mup_LastStationT->clear(); mup_OneStationT->clear();
      mup_LastStationAngT->clear(); mup_OneStationAngT->clear(); mup_2DCompatibilityT->clear();
      mup_NMuonHits->clear(); mup_NMuonStations->clear(); mup_NTrackerLayers->clear(); mup_NPixelLayers->clear(); mup_relIso->clear();

      massPiPiPi->clear(); pfMass->clear();
      pionPF_Eta->clear(); pionPF_Phi->clear();
      pionPF_Pt->clear();  pi1_maxDeltaR->clear(); pi2_maxDeltaR->clear(); pi1_Pt->clear();
      pi2_Pt->clear(); pi1_numberOfGammas->clear();
      pi2_numberOfGammas->clear();
   }



   int const Bfinder::getMuCat(reco::Muon const& muon) const{
     int muCat = 0;
     if (muon.isGlobalMuon()) {
       if (muon.isTrackerMuon()) muCat = 1;
       else muCat = 10;
     }
     else if (muon.isTrackerMuon()) {
       if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)) {
         if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)) muCat = 6;
         else if (muon::isGoodMuon(muon, muon::TMOneStationTight)) muCat = 5;
         else if (muon::isGoodMuon(muon, muon::TMOneStationLoose)) muCat = 4;
         else muCat = 3;
       }
       else muCat = 2;
     }
     else if (muon.isStandAloneMuon()) muCat = 7;
     else if (muon.isCaloMuon()) muCat = 8;
     else muCat = 9;

     if ( !(muon::isGoodMuon(muon, muon::TMOneStationLoose)) && muon::isGoodMuon(muon, muon::TMOneStationTight) )
       std::cout << "inconsistent muon cat 1" << std::endl;
     if ( !(muon::isGoodMuon(muon, muon::TMOneStationTight)) && muon::isGoodMuon(muon, muon::TMLastStationAngTight) )
       std::cout << "inconsistent muon cat 2" << std::endl;

     return muCat;
   }



   // ------------ method called once each job just before starting event loop  ------------

   void Bfinder::beginJob()
   {

     std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

     edm::Service<TFileService> fs;
     tree_ = fs->make<TTree>("ntuple","B ntuple");

     tree_->Branch("nB"                , &nB       , "nB/i"    );
     tree_->Branch("nMu"               , &nMu      , "nMu/i"   );
     tree_->Branch("nVtx"              , &nVtx     , "nVtx/i"  );

     tree_->Branch("psiMass"            , &psiMass               );

     tree_->Branch("B_J_mass"          , &B_J_mass             );
     tree_->Branch("B_J_px"            , &B_J_px               );
     tree_->Branch("B_J_py"            , &B_J_py               );
     tree_->Branch("B_J_pz"            , &B_J_pz               );
     tree_->Branch("B_J_DecayVtxX"     , &B_J_DecayVtxX        );
     tree_->Branch("B_J_DecayVtxY"     , &B_J_DecayVtxY        );
     tree_->Branch("B_J_DecayVtxZ"     , &B_J_DecayVtxZ        );
     tree_->Branch("B_J_DecayVtxXE"    , &B_J_DecayVtxXE       );
     tree_->Branch("B_J_DecayVtxYE"    , &B_J_DecayVtxYE       );
     tree_->Branch("B_J_DecayVtxZE"    , &B_J_DecayVtxZE       );

     tree_->Branch("B_mu_px1_cjp"      , &B_mu_px1_cjp         );
     tree_->Branch("B_mu_py1_cjp"      , &B_mu_py1_cjp         );
     tree_->Branch("B_mu_pz1_cjp"      , &B_mu_pz1_cjp         );

     tree_->Branch("B_mu_px2_cjp"      , &B_mu_px2_cjp         );
     tree_->Branch("B_mu_py2_cjp"      , &B_mu_py2_cjp         );
     tree_->Branch("B_mu_pz2_cjp"      , &B_mu_pz2_cjp         );

     tree_->Branch("B_J_Prob"          , &B_J_Prob             );

     tree_->Branch("kaonP_px_0c"        , &kaonP_px_0c           );
     tree_->Branch("kaonP_py_0c"        , &kaonP_py_0c           );
     tree_->Branch("kaonP_pz_0c"        , &kaonP_pz_0c           );
     tree_->Branch("kaonP_pt_0c"        , &kaonP_pt_0c           );

     tree_->Branch("kaonP_track_normchi2"   , &kaonP_track_normchi2      );
     tree_->Branch("kaonP_Hits"        , &kaonP_Hits           );
     tree_->Branch("kaonP_PHits"       , &kaonP_PHits          );
     tree_->Branch("kaonP_NTrackerLayers"       , &kaonP_NTrackerLayers          );
     tree_->Branch("kaonP_NPixelLayers"       , &kaonP_NPixelLayers          );

     tree_->Branch("kaonM_px_0c"        , &kaonM_px_0c           );
     tree_->Branch("kaonM_py_0c"        , &kaonM_py_0c           );
     tree_->Branch("kaonM_pz_0c"        , &kaonM_pz_0c           );
     tree_->Branch("kaonM_pt_0c"        , &kaonM_pt_0c           );

     tree_->Branch("kaonM_track_normchi2"   , &kaonM_track_normchi2      );
     tree_->Branch("kaonM_Hits"        , &kaonM_Hits           );
     tree_->Branch("kaonM_PHits"       , &kaonM_PHits          );
     tree_->Branch("kaonM_NTrackerLayers"       , &kaonM_NTrackerLayers          );
     tree_->Branch("kaonM_NPixelLayers" , &kaonM_NPixelLayers );

     tree_->Branch("B_Phi_Prob" , &B_Phi_Prob );
     tree_->Branch("maxDelta" , &maxDelta);



     tree_->Branch("PV_bestBang_X"     , &PV_bestBang_X        );
     tree_->Branch("PV_bestBang_Y"     , &PV_bestBang_Y        );
     tree_->Branch("PV_bestBang_Z"     , &PV_bestBang_Z        );
     tree_->Branch("PV_bestBang_XE"    , &PV_bestBang_XE       );
     tree_->Branch("PV_bestBang_YE"    , &PV_bestBang_YE       );
     tree_->Branch("PV_bestBang_ZE"    , &PV_bestBang_ZE       );
     tree_->Branch("PV_bestBang_XYE"   , &PV_bestBang_XYE      );
     tree_->Branch("PV_bestBang_XZE"   , &PV_bestBang_XZE      );
     tree_->Branch("PV_bestBang_YZE"   , &PV_bestBang_YZE      );
     tree_->Branch("PV_bestBang_CL"    , &PV_bestBang_CL       );

     tree_->Branch("PV_bestBang_RF_X"  , &PV_bestBang_RF_X     );
     tree_->Branch("PV_bestBang_RF_Y"  , &PV_bestBang_RF_Y     );
     tree_->Branch("PV_bestBang_RF_Z"  , &PV_bestBang_RF_Z     );
     tree_->Branch("PV_bestBang_RF_XE" , &PV_bestBang_RF_XE    );
     tree_->Branch("PV_bestBang_RF_YE" , &PV_bestBang_RF_YE    );
     tree_->Branch("PV_bestBang_RF_ZE" , &PV_bestBang_RF_ZE    );
     tree_->Branch("PV_bestBang_RF_XYE", &PV_bestBang_RF_XYE   );
     tree_->Branch("PV_bestBang_RF_XZE", &PV_bestBang_RF_XZE   );
     tree_->Branch("PV_bestBang_RF_YZE", &PV_bestBang_RF_YZE   );
     tree_->Branch("PV_bestBang_RF_CL" , &PV_bestBang_RF_CL    );
     tree_->Branch("PV_bestBang_RF_NTrkDif",&PV_bestBang_RF_NTrkDif  );

     tree_->Branch("run"               , &run                  , "run/I"                   );
     tree_->Branch("event"             , &event                , "event/I"                 );

     tree_->Branch("nTrgL"             , &nTrgL                , "nTrgL/I"                 );
     tree_->Branch("nTrgL1L"           , &nTrgL1L              , "nTrgL1L/I"               );
     tree_->Branch("triggersL"         , &triggersL            , "triggersL[nTrgL]/C"      );
   //   tree_->Branch("triggersL1L"       , &triggersL1L          , "triggersL1L[nTrgL1L]/C"  );
     tree_->Branch("triggersMuPL"      , &triggersMuPL         );
     tree_->Branch("triggersMuML"      , &triggersMuML         );
     tree_->Branch("triggersL1L2_MuPL" , &triggersL1L2_MuPL    );
     tree_->Branch("triggersL1L2_MuML" , &triggersL1L2_MuML    );

   //////------------////////
     tree_->Branch("mum_normChi2"      , &mum_normChi2         );
     tree_->Branch("mumCat"            , &mumCat               );
     tree_->Branch("mumAngT"           , &mumAngT              );
     tree_->Branch("mumNHits"          , &mumNHits             );
     tree_->Branch("mumNPHits"         , &mumNPHits            );
     tree_->Branch("mum_dxy_Bsdecay"            , &mum_dxy_Bsdecay               );
     tree_->Branch("mum_dz_Bsdecay"             , &mum_dz_Bsdecay                );

     tree_->Branch("mum_isGlobalMuon"               , &mum_isGlobalMuon     );
     tree_->Branch("mum_isTrackerMuon"              , &mum_isTrackerMuon    );
     tree_->Branch("mum_isTight"                    , &mum_isTight          );
     tree_->Branch("mum_LastStationOptimizedLowPtT" , &mum_LastStationOptimizedLowPtT  );
     tree_->Branch("mum_LastStationT"               , &mum_LastStationT          );
     tree_->Branch("mum_OneStationT"                , &mum_OneStationT          );
     tree_->Branch("mum_LastStationAngT"            , &mum_LastStationAngT          );
     tree_->Branch("mum_OneStationAngT"             , &mum_OneStationAngT          );
     tree_->Branch("mum_2DCompatibilityT"           , &mum_2DCompatibilityT          );

     tree_->Branch("mum_NMuonHits"     , &mum_NMuonHits        );
     tree_->Branch("mum_NMuonStations" , &mum_NMuonStations    );
     tree_->Branch("mum_NTrackerLayers", &mum_NTrackerLayers   );
     tree_->Branch("mum_NPixelLayers"  , &mum_NPixelLayers     );
     tree_->Branch("mum_relIso"        , &mum_relIso           );


   //////------------////////
     tree_->Branch("mup_normChi2"      , &mup_normChi2         );
     tree_->Branch("mupCat"            , &mupCat               );
     tree_->Branch("mupAngT"           , &mupAngT              );
     tree_->Branch("mupNHits"          , &mupNHits             );
     tree_->Branch("mupNPHits"         , &mupNPHits            );
     tree_->Branch("mup_dxy_Bsdecay"            , &mup_dxy_Bsdecay               );
     tree_->Branch("mup_dz_Bsdecay"             , &mup_dz_Bsdecay                );

     tree_->Branch("mup_isGlobalMuon"    , &mup_isGlobalMuon     );
     tree_->Branch("mup_isTrackerMuon"   , &mup_isTrackerMuon    );
     tree_->Branch("mup_isTight"         , &mup_isTight          );
     tree_->Branch("mup_LastStationOptimizedLowPtT" , &mup_LastStationOptimizedLowPtT  );
     tree_->Branch("mup_LastStationT"               , &mup_LastStationT          );
     tree_->Branch("mup_OneStationT"                , &mup_OneStationT          );
     tree_->Branch("mup_LastStationAngT"            , &mup_LastStationAngT          );
     tree_->Branch("mup_OneStationAngT"             , &mup_OneStationAngT          );
     tree_->Branch("mup_2DCompatibilityT"           , &mup_2DCompatibilityT          );


     tree_->Branch("mup_NMuonHits"     , &mup_NMuonHits        );
     tree_->Branch("mup_NMuonStations" , &mup_NMuonStations    );
     tree_->Branch("mup_NTrackerLayers", &mup_NTrackerLayers   );
     tree_->Branch("mup_NPixelLayers"  , &mup_NPixelLayers     );
     tree_->Branch("mup_relIso"        , &mup_relIso           );

     tree_->Branch("massPiPiPi"               , &massPiPiPi          );
     tree_->Branch("pfMass"               , &pfMass          );

     tree_->Branch("pionPF_Eta"               , &pionPF_Eta          );
     tree_->Branch("pionPF_Phi"               , &pionPF_Phi          );
     tree_->Branch("pionPF_Pt"               , &pionPF_Pt          );
     tree_->Branch("pi1_maxDeltaR"               , &pi1_maxDeltaR          );
     tree_->Branch("pi2_maxDeltaR"               , &pi2_maxDeltaR          );
     tree_->Branch("pi1_Pt"               , &pi1_Pt          );
     tree_->Branch("pi2_Pt"               , &pi2_Pt          );
     tree_->Branch("pi1_numberOfGammas"               , &pi1_numberOfGammas          );
     tree_->Branch("pi2_numberOfGammas"               , &pi2_numberOfGammas          );

   }


   // ------------ method called once each job just after ending the event loop  ------------
   void Bfinder::endJob() {
     tree_->GetDirectory()->cd();
     tree_->Write();
   }

   //define this as a plug-in
   DEFINE_FWK_MODULE(Bfinder);
