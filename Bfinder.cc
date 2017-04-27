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

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "TLorentzVector.h"
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

   B_mass_cjp(0),
   B_px_cjp(0),           B_py_cjp(0),          B_pz_cjp(0),
   B_DecayVtxX(0),    B_DecayVtxY(0),   B_DecayVtxZ(0),
   B_DecayVtxXE(0),   B_DecayVtxYE(0),  B_DecayVtxZE(0),


   B_J_mass(0),       B_J_px(0),            B_J_py(0),        B_J_pz(0),
   B_J_DecayVtxX(0),  B_J_DecayVtxY(0),     B_J_DecayVtxZ(0),
   B_J_DecayVtxXE(0), B_J_DecayVtxYE(0),    B_J_DecayVtxZE(0),

   B_pion_px_cjp(0),   B_pion_py_cjp(0),  B_pion_pz_cjp(0),
   B_pion_charge(0),   pion_track_normchi2(0),     pion_Hits(0),  pion_PHits(0),
   pion_dxy_Bcdecay(0), pion_dz_Bcdecay(0), pion_NTrackerLayers(0),  pion_NPixelLayers(0), pi_Bcdecay_weight(0),

   B_mu_px1_cjp(0),   B_mu_py1_cjp(0),  B_mu_pz1_cjp(0),
   B_mu_px2_cjp(0),   B_mu_py2_cjp(0),  B_mu_pz2_cjp(0),

   B_Prob(0), B_J_Prob(0),

   PV_bestBang_X(0),      PV_bestBang_Y(0),     PV_bestBang_Z(0),
   PV_bestBang_XE(0),     PV_bestBang_YE(0),    PV_bestBang_ZE(0),
   PV_bestBang_XYE(0),    PV_bestBang_XZE(0),   PV_bestBang_YZE(0),
   PV_bestBang_CL(0),

   PV_bestBang_RF_X(0),   PV_bestBang_RF_Y(0),  PV_bestBang_RF_Z(0),
   PV_bestBang_RF_XE(0),  PV_bestBang_RF_YE(0), PV_bestBang_RF_ZE(0),
   PV_bestBang_RF_XYE(0), PV_bestBang_RF_XZE(0),PV_bestBang_RF_YZE(0),
   PV_bestBang_RF_CL(0),  PV_bestBang_RF_NTrkDif(0),

   mum_normChi2(0),  mumdxy(0),    mumdz(0)    , mumCat(0) , mumAngT(0)    , mumNHits(0)   , mumNPHits(0),
   mum_isGlobalMuon(0), mum_isTight(0), mum_NMuonHits(0), mum_NMuonStations(0), mum_NTrackerLayers(0), mum_NPixelLayers(0), mum_relIso(0),

   mup_normChi2(0),  mupdxy(0),    mupdz(0)    , mupCat(0) , mupAngT(0)    , mupNHits(0)   , mupNPHits(0),
   mup_isGlobalMuon(0), mup_isTight(0), mup_NMuonHits(0), mup_NMuonStations(0), mup_NTrackerLayers(0), mup_NPixelLayers(0), mup_relIso(0),
   phoEta(0),
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
  iEvent.getByLabel("photons", photonHandle);
//  iEvent.getByLabel(photonCollectionProducer_, photonCollection_ , photonHandle);
//  const reco::PhotonCollection photonCollection = *(photonHandle.product());


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
    std::string::size_type trigger1 = alltnames.find("HLT_DoubleMu4_Jpsi_Displaced",0);
//    std::string::size_type trigger1 = alltnames.find("HLT_DoubleMu3p5_LowMass_Displaced",0);

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
      ParticleMass PM_PDG_JPSI_MASS = PDG_JPSI_MASS;
      ParticleMass PM_PDG_PION_MASS = PDG_PION_MASS;
      float PM_muon_sigma = PM_PDG_MUON_MASS*1.e-6;
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
         triggersMuPL ->push_back(triggersMuP);
         triggersMuML ->push_back(triggersMuM);


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


	   pat::GenericParticle patTrack1;
// 	   pat::GenericParticle patTrack2;

	   for(vector<pat::GenericParticle>::const_iterator iTrack1 = thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end(); ++iTrack1 )
	     {

// 	       int ngenT1 = 0;//PdgIDatTruthLevel(iTrack1->track(), genParticles, PId1);
	       patTrack1 = *iTrack1;
	       if(iTrack1->pt()<2.0) continue;
           if (fabs(iTrack1->eta())>2.5) continue;

	       if(!(patTrack1.track()->quality(reco::TrackBase::highPurity))) continue;

               bool matchflag = false;
               const reco::CandidatePtrVector & mu1P_overlaps = patTrack1.overlaps(muonTypeForPAT);
               if ( mu1P_overlaps.size() > 0 ) //std::cout << "patTrack1 overlaps with a muon." << endl;
               for (size_t i = 0; i < mu1P_overlaps.size(); ++i) {
                 const pat::Muon *mu = dynamic_cast<const pat::Muon *>(&*mu1P_overlaps[i]);
                 if (mu) {
                   // check here that muon match isn't the same as a muon used in the reco...
                   if (mu==patMuonP || mu==patMuonM)  {
                       //std::cout << "match between patTrack1 and patMuonP/patMuonM" << endl;
                       matchflag=true;
                   }
                 }
               }

               if(matchflag) continue;

               TransientTrack pion1TT(patTrack1.track(), &(*bFieldHandle) );

               TLorentzVector p4k1_0c, p4mumuk_0c;
               p4k1_0c.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(), PDG_PION_MASS);

               p4mumuk_0c = p4mu1_0c + p4mu2_0c + p4k1_0c;
               if(p4mumuk_0c.M() > 7.0) continue;
               if(p4mumuk_0c.M() < 5.5) continue;

               //Now we are ready to combine!

               std::vector<RefCountedKinematicParticle> B_candidate_init;

               B_candidate_init.push_back(pFactory.particle(muon1TT, PM_PDG_MUON_MASS, chi,ndf, PM_muon_sigma));
               B_candidate_init.push_back(pFactory.particle(muon2TT, PM_PDG_MUON_MASS, chi,ndf, PM_muon_sigma));
               B_candidate_init.push_back(pFactory.particle(pion1TT, PM_PDG_PION_MASS, chi,ndf, PM_muon_sigma));
               RefCountedKinematicTree xbVFT, vertexFitTree;

               std::vector<RefCountedKinematicParticle> B_candidate = B_candidate_init;
               KinematicParticleVertexFitter pFitter; //KinematicParticleVertexFitter
               xbVFT = pFitter.fit(B_candidate);

               if (!xbVFT->isValid()) continue;
               xbVFT->movePointerToTheTop();


    ///////////////////////////////  BS VERTEX FIT /////////////////////////////////////////
               //J/Psi mass constraint

               B_candidate = B_candidate_init;
               MultiTrackKinematicConstraint *ConstraintJpsiMass = new TwoTrackMassKinematicConstraint(PM_PDG_JPSI_MASS);
               KinematicConstrainedVertexFitter kcvFitter; //KinematicParticleVertexFitter
               vertexFitTree = kcvFitter.fit(B_candidate, ConstraintJpsiMass);

               if (!vertexFitTree->isValid()) continue;

               vertexFitTree->movePointerToTheTop();
               RefCountedKinematicParticle bCandCjp = vertexFitTree->currentParticle();
               RefCountedKinematicVertex bDecayVertexCjp = vertexFitTree->currentDecayVertex();
               if (!bDecayVertexCjp->vertexIsValid())  continue;

               double B_mass_cjp_tmp = bCandCjp->currentState().mass();
               if(B_mass_cjp_tmp < 6.0) continue;
               if(B_mass_cjp_tmp > 6.6) continue;
               //
               if(bDecayVertexCjp->chiSquared()<0) continue;
               double B_Prob_tmp   = TMath::Prob(bDecayVertexCjp->chiSquared(), (int) bDecayVertexCjp->degreesOfFreedom());
               if(B_Prob_tmp < 0.01) continue;

               // get children from final B fit

               vertexFitTree->movePointerToTheFirstChild();
               RefCountedKinematicParticle mu1CandMC    = vertexFitTree->currentParticle();
               vertexFitTree->movePointerToTheNextChild();
               RefCountedKinematicParticle mu2CandMC    = vertexFitTree->currentParticle();
               vertexFitTree->movePointerToTheNextChild();
               RefCountedKinematicParticle T1CandMC     = vertexFitTree->currentParticle();


               ROOT::Math::XYZPoint bDecayPoint( (*bDecayVertexCjp).position().x(), (*bDecayVertexCjp).position().y(), (*bDecayVertexCjp).position().z() );
 	       GlobalPoint BcGP = GlobalPoint( (*bDecayVertexCjp).position().x(), (*bDecayVertexCjp).position().y(), (*bDecayVertexCjp).position().z() );


        // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE BS TRACKS FROM ITS FIT
      // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

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

                    Double_t dx = (*bDecayVertexCjp).position().x() - PVtxBeSp.x();
                    Double_t dy = (*bDecayVertexCjp).position().y() - PVtxBeSp.y();
                    Double_t dz = (*bDecayVertexCjp).position().z() - PVtxBeSp.z();
                    Double_t cosAlphaXYZ = ( bCandCjp->currentState().globalMomentum().x() * dx + bCandCjp->currentState().globalMomentum().y()*dy + bCandCjp->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* bCandCjp->currentState().globalMomentum().mag() );

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

     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     /////////////////////////////////// // try refitting the primary without the tracks in the B reco candidate

               // first get tracks from the original primary
               vector<reco::TransientTrack> vertexTracks;

               for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestPV_Bang.tracks_begin();
                 iTrack != bestPV_Bang.tracks_end(); ++iTrack)
              {
                    // compare primary tracks to check for matches with B cand
                    TrackRef trackRef = iTrack->castTo<TrackRef>();

                    // the 4 tracks in the B cand are  patTrack1 glbTrackP glbTrackM
                    if (  !(   (patTrack1.track()==trackRef)   ||
                               (trkTrackP==trackRef)           ||
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
                  //bctauRf_temp = Myctau(bCandCjp, bDecayVertexCjp, PVRfP, PVRfE, mb, bctau2DRf_temp, bctauRfE_temp, bctau2DRfE_temp);

                  //set bestVtxRf as new best vertex to fill variables for ntuple
                  bestVtxRf = reco::Vertex(v);
                }
               }



  for( std::vector<reco::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho)
    {
        reco::Photon localPho = reco::Photon(*iPho);
        phoEta->push_back( localPho.eta() );

        ///~~~fit 3 tracks from Bc together~~~///
        	vector<reco::TransientTrack> BcTracks;
        	BcTracks.push_back(muon1TT);
        	BcTracks.push_back(muon2TT);
        	BcTracks.push_back(pion1TT);

                AdaptiveVertexFitter BcVertexFitter;
                TransientVertex BcTV = BcVertexFitter.vertex(BcTracks, BcGP);

                if ( BcTV.isValid() )
        	{

                   BcVtx = reco::Vertex(BcTV);
        	   pi_Bcdecay_weight->push_back(BcVtx.trackWeight(patTrack1.track()));
                }
        	else
        	{
        	   pi_Bcdecay_weight->push_back(-1);
        	}


          PV_bestBang_RF_NTrkDif->push_back( bestPV_Bang.tracksSize() - vertexTracks.size() );

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

               //cout << "PV bestPV_Bang: " <<bestPV_Bang.x()<< " "<<bestPV_Bang.y()<<" "<<bestPV_Bang.z()<< endl;

               GlobalVector mu1CandMC_p = mu1CandMC->currentState().kinematicParameters().momentum();
               GlobalVector mu2CandMC_p = mu2CandMC->currentState().kinematicParameters().momentum();

               if(iMuon1->charge() == 1 ) mupCategory = getMuCat( *iMuon1 );
               if(iMuon1->charge() == -1) mumCategory = getMuCat( *iMuon1 );
               if(iMuon2->charge() == 1 ) mupCategory = getMuCat( *iMuon2 );
               if(iMuon2->charge() == -1) mumCategory = getMuCat( *iMuon2 );

               const reco::Muon *recoMuonM = patMuonM;
               const reco::Muon *recoMuonP = patMuonP;

               // fill candidate variables now

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


    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// HERE INDEX IS [nB] /////////////////////////
               B_mass_cjp           ->push_back(B_mass_cjp_tmp);
               B_px_cjp             ->push_back(bCandCjp->currentState().globalMomentum().x());
               B_py_cjp             ->push_back(bCandCjp->currentState().globalMomentum().y());
               B_pz_cjp             ->push_back(bCandCjp->currentState().globalMomentum().z());
               B_DecayVtxX      ->push_back(bDecayVertexCjp->position().x());
               B_DecayVtxY      ->push_back(bDecayVertexCjp->position().y());
               B_DecayVtxZ      ->push_back(bDecayVertexCjp->position().z());
               B_DecayVtxXE     ->push_back(bDecayVertexCjp->error().cxx());
               B_DecayVtxYE     ->push_back(bDecayVertexCjp->error().cyy());
               B_DecayVtxZE     ->push_back(bDecayVertexCjp->error().czz());
//                B_DecayVtxXYE->push_back(bDecayVertexCjp->error().cyx());
//                B_DecayVtxXZE->push_back(bDecayVertexCjp->error().czx());
//                B_DecayVtxYZE->push_back(bDecayVertexCjp->error().czy());


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

               B_pion_px_cjp       ->push_back(T1CandMC->currentState().globalMomentum().x());
               B_pion_py_cjp       ->push_back(T1CandMC->currentState().globalMomentum().y());
               B_pion_pz_cjp       ->push_back(T1CandMC->currentState().globalMomentum().z());
               B_pion_charge   ->push_back(iTrack1->charge());
               pion_track_normchi2  ->push_back(patTrack1.track()->normalizedChi2());
               pion_Hits       ->push_back(patTrack1.track()->numberOfValidHits() );
               pion_PHits      ->push_back(patTrack1.track()->hitPattern().numberOfValidPixelHits() );
               pion_dxy_Bcdecay	->push_back(patTrack1.track()->dxy(bDecayPoint) );
               pion_dz_Bcdecay		->push_back(patTrack1.track()->dz(bDecayPoint	) );
	       pion_NTrackerLayers->push_back ( patTrack1.track()->hitPattern().trackerLayersWithMeasurement() );
	       pion_NPixelLayers->push_back ( patTrack1.track()->hitPattern().pixelLayersWithMeasurement() );


               B_mu_px1_cjp ->push_back(mu1CandMC_p.x());
               B_mu_py1_cjp ->push_back(mu1CandMC_p.y());
               B_mu_pz1_cjp ->push_back(mu1CandMC_p.z());
//                B_J_charge1->push_back(mНУ А как же всем известный АМХ 40 ?)))u1CandMC->currentState().particleCharge()); // noneed

               B_mu_px2_cjp ->push_back(mu2CandMC_p.x());
               B_mu_py2_cjp ->push_back(mu2CandMC_p.y());
               B_mu_pz2_cjp ->push_back(mu2CandMC_p.z());
//                B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());


               B_Prob    ->push_back(B_Prob_tmp);
               B_J_Prob  ->push_back(JP_Prob_tmp);

//------------------//
               mumCat->push_back( mumCategory );
	       mum_isGlobalMuon->push_back ( recoMuonM->isGlobalMuon() );
               mumAngT->push_back( muon::isGoodMuon(*recoMuonM, muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
     	       mum_isTight->push_back ( muon::isTightMuon(*recoMuonM, BcVtx ) ); // Might not be correct
//     	       mum_isTight->push_back (0);

               mumdxy->push_back( recoMuonM->muonBestTrack()->dxy(bDecayPoint) );// el dxy del Muon negatico respetcto del PV con BSc (el de mayor pt)
               mumdz->push_back( recoMuonM->muonBestTrack()->dz(bDecayPoint) );

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

		mum_NTrackerLayers->push_back ( recoMuonM->muonBestTrack()->hitPattern().trackerLayersWithMeasurement() );
		mum_NPixelLayers->push_back ( recoMuonM->muonBestTrack()->hitPattern().pixelLayersWithMeasurement() );
                mumNPHits->push_back( recoMuonM->muonBestTrack()->hitPattern().numberOfValidPixelHits() );
                mumNHits->push_back( recoMuonM->muonBestTrack()->numberOfValidHits() );

		float mum_chIso  = patMuonM->chargedHadronIso();
		float mum_nhIso  = patMuonM->neutralHadronIso();
		float mum_phIso  = patMuonM->photonIso();
		float mum_puIso  = patMuonM->puChargedHadronIso();
		mum_relIso->push_back( (mum_chIso + std::max(0.0, mum_nhIso + mum_phIso - 0.5 * mum_puIso)) / patMuonM->pt() );


//------------------//
               mupCat->push_back( mupCategory );
	       mup_isGlobalMuon->push_back ( recoMuonP->isGlobalMuon() );
               mupAngT->push_back( muon::isGoodMuon(*recoMuonP, muon::TMOneStationTight) ); // este es para poner la condicion si es o no softmuon
     	       mup_isTight->push_back ( muon::isTightMuon(*recoMuonP, BcVtx) );  // Might not be correct
//     	       mup_isTight->push_back (0);

               mupdxy->push_back( recoMuonP->muonBestTrack()->dxy(bDecayPoint) );// el dxy del Muon negatico respetcto del PV con BSc (el de mayor pt)
               mupdz->push_back( recoMuonP->muonBestTrack()->dz(bDecayPoint) );

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

		mup_NTrackerLayers->push_back ( recoMuonP->muonBestTrack()->hitPattern().trackerLayersWithMeasurement() );
		mup_NPixelLayers->push_back ( recoMuonP->muonBestTrack()->hitPattern().pixelLayersWithMeasurement() );
                mupNPHits->push_back( recoMuonP->muonBestTrack()->hitPattern().numberOfValidPixelHits() );
                mupNHits->push_back( recoMuonP->muonBestTrack()->numberOfValidHits() );

		float mup_chIso  = patMuonP->chargedHadronIso();
		float mup_nhIso  = patMuonP->neutralHadronIso();
		float mup_phIso  = patMuonP->photonIso();
		float mup_puIso  = patMuonP->puChargedHadronIso();
		mup_relIso->push_back( (mup_chIso + std::max(0.0, mup_nhIso + mup_phIso - 0.5 * mup_puIso)) / patMuonP->pt() );

               nB++;

               /////////////////////////////////////////////////

               muonParticles.clear();
               B_candidate_init.clear();
               B_candidate.clear();
	     }
	}
    }
}

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

   B_mass_cjp->clear();
   B_px_cjp->clear();           B_py_cjp->clear();          B_pz_cjp->clear();
   B_DecayVtxX->clear();    B_DecayVtxY->clear();   B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();   B_DecayVtxYE->clear();  B_DecayVtxZE->clear();

   B_J_mass->clear();       B_J_px->clear();            B_J_py->clear();        B_J_pz->clear();
   B_J_DecayVtxX->clear();  B_J_DecayVtxY->clear();     B_J_DecayVtxZ->clear();
   B_J_DecayVtxXE->clear(); B_J_DecayVtxYE->clear();    B_J_DecayVtxZE->clear();

   B_pion_px_cjp->clear();     B_pion_py_cjp->clear();    B_pion_pz_cjp->clear();
   B_pion_charge->clear(); pion_track_normchi2->clear();   pion_Hits->clear();    pion_PHits->clear();
   pion_dxy_Bcdecay->clear();  pion_dz_Bcdecay->clear(); pion_NTrackerLayers->clear();  pion_NPixelLayers->clear();
   pi_Bcdecay_weight->clear();

   B_mu_px1_cjp->clear();   B_mu_py1_cjp->clear();  B_mu_pz1_cjp->clear();
   B_mu_px2_cjp->clear();   B_mu_py2_cjp->clear();  B_mu_pz2_cjp->clear();

   B_Prob->clear(); B_J_Prob->clear();

   PV_bestBang_X->clear();      PV_bestBang_Y->clear();     PV_bestBang_Z->clear();
   PV_bestBang_XE->clear();     PV_bestBang_YE->clear();    PV_bestBang_ZE->clear();
   PV_bestBang_XYE->clear();    PV_bestBang_XZE->clear();   PV_bestBang_YZE->clear();
   PV_bestBang_CL->clear();

   PV_bestBang_RF_X->clear();   PV_bestBang_RF_Y->clear();  PV_bestBang_RF_Z->clear();
   PV_bestBang_RF_XE->clear();  PV_bestBang_RF_YE->clear(); PV_bestBang_RF_ZE->clear();
   PV_bestBang_RF_XYE->clear(); PV_bestBang_RF_XZE->clear();PV_bestBang_RF_YZE->clear();
   PV_bestBang_RF_CL->clear();  PV_bestBang_RF_NTrkDif->clear();

   mum_normChi2->clear();  mumdxy->clear();    mumdz->clear(); mumCat->clear();    mumAngT->clear();   mumNHits->clear();  mumNPHits->clear();
   mum_isGlobalMuon->clear(); mum_isTight->clear(); mum_NMuonHits->clear(); mum_NMuonStations->clear(); mum_NTrackerLayers->clear(); mum_NPixelLayers->clear(); mum_relIso->clear();

   mup_normChi2->clear();   mupdxy->clear();    mupdz->clear(); mupCat->clear();    mupAngT->clear();   mupNHits->clear();  mupNPHits->clear();
   mup_isGlobalMuon->clear(); mup_isTight->clear(); mup_NMuonHits->clear(); mup_NMuonStations->clear(); mup_NTrackerLayers->clear(); mup_NPixelLayers->clear(); mup_relIso->clear();
   phoEta->clear();

   // }}}
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

  tree_->Branch("B_mass_cjp"            , &B_mass_cjp               );
  tree_->Branch("B_px_cjp"              , &B_px_cjp                 );
  tree_->Branch("B_py_cjp"              , &B_py_cjp                 );
  tree_->Branch("B_pz_cjp"              , &B_pz_cjp                 );
  tree_->Branch("B_DecayVtxX"       , &B_DecayVtxX          );
  tree_->Branch("B_DecayVtxY"       , &B_DecayVtxY          );
  tree_->Branch("B_DecayVtxZ"       , &B_DecayVtxZ          );
  tree_->Branch("B_DecayVtxXE"      , &B_DecayVtxXE         );
  tree_->Branch("B_DecayVtxYE"      , &B_DecayVtxYE         );
  tree_->Branch("B_DecayVtxZE"      , &B_DecayVtxZE         );


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

  tree_->Branch("B_pion_px_cjp"        , &B_pion_px_cjp           );
  tree_->Branch("B_pion_py_cjp"        , &B_pion_py_cjp           );
  tree_->Branch("B_pion_pz_cjp"        , &B_pion_pz_cjp           );
  tree_->Branch("B_pion_charge"        , &B_pion_charge       	  );
  tree_->Branch("pion_track_normchi2"  , &pion_track_normchi2     );
  tree_->Branch("pion_Hits"            , &pion_Hits               );
  tree_->Branch("pion_PHits"           , &pion_PHits              );
  tree_->Branch("pion_dxy_Bcdecay"     , &pion_dxy_Bcdecay        );
  tree_->Branch("pion_dz_Bcdecay"      , &pion_dz_Bcdecay         );
  tree_->Branch("pion_NTrackerLayers"  , &pion_NTrackerLayers     );
  tree_->Branch("pion_NPixelLayers"    , &pion_NPixelLayers       );
  tree_->Branch("pi_Bcdecay_weight"    , &pi_Bcdecay_weight       );

  tree_->Branch("B_mu_px1_cjp"      , &B_mu_px1_cjp         );
  tree_->Branch("B_mu_py1_cjp"      , &B_mu_py1_cjp         );
  tree_->Branch("B_mu_pz1_cjp"      , &B_mu_pz1_cjp         );

  tree_->Branch("B_mu_px2_cjp"      , &B_mu_px2_cjp         );
  tree_->Branch("B_mu_py2_cjp"      , &B_mu_py2_cjp         );
  tree_->Branch("B_mu_pz2_cjp"      , &B_mu_pz2_cjp         );

  tree_->Branch("B_Prob"            , &B_Prob               );
  tree_->Branch("B_J_Prob"          , &B_J_Prob             );

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
  tree_->Branch("mumdxy"            , &mumdxy               );
  tree_->Branch("mumdz"             , &mumdz                );

  tree_->Branch("mum_isGlobalMuon"  , &mum_isGlobalMuon     );
  tree_->Branch("mum_isTight"       , &mum_isTight          );
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
  tree_->Branch("mupdxy"            , &mupdxy               );
  tree_->Branch("mupdz"             , &mupdz                );

  tree_->Branch("mup_isGlobalMuon"  , &mup_isGlobalMuon     );
  tree_->Branch("mup_isTight"       , &mup_isTight          );
  tree_->Branch("mup_NMuonHits"     , &mup_NMuonHits        );
  tree_->Branch("mup_NMuonStations" , &mup_NMuonStations    );
  tree_->Branch("mup_NTrackerLayers", &mup_NTrackerLayers   );
  tree_->Branch("mup_NPixelLayers"  , &mup_NPixelLayers     );
  tree_->Branch("mup_relIso"        , &mup_relIso           );
  tree_->Branch("phoEta"            , &phoEta           );

}


// ------------ method called once each job just after ending the event loop  ------------
void Bfinder::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bfinder);
