import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")


from PhysicsTools.PatAlgos.patTemplate_cfg import *


process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1), SkipEvent = cms.untracked.vstring('ProductNotFound') )
#process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(1000), SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(

'/store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20000/16E502D6-7068-E211-B692-001A645C1EC4.root'
#'/store/data/Run2012C/MuOniaParked/AOD/22Jan2013-v1/20000/00839611-2576-E211-B384-E41F13181568.root'

#'/store/data/Run2012B/MuOnia/AOD/22Jan2013-v1/20000/0215100D-C784-E211-A2DD-002618943972.root'
#'/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/000F9808-D483-E211-975C-003048FFD7BE.root'
#'/store/data/Run2012D/MuOnia/AOD/22Jan2013-v1/10000/004C60F9-188F-E211-B43F-001E6849D384.root'
	)
)

#from myAnalyzers.JPsiKsPAT.RecoInput2_cfi import *

process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START311_V2::All')
# process.GlobalTag.globaltag = cms.string('FT_R_53_V18::All') ## their
process.GlobalTag.globaltag = cms.string( 'FT53_V21A_AN6::All') ## our for 2012ABC
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.GlobalTag.globaltag = 'START44_V9B::All'

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,					    #	      patAODTrackCands
	label='TrackCands',		      # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
	tracks=cms.InputTag('generalTracks'), # input track collection
	particleType='pi+',		      # particle type (for assigning a mass)
	preselection='pt > 0.4',	      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
	selection='pt > 0.4',		      # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
	isolation={},			      # Isolations to use ('source':deltaR; set to {} for None)
	isoDeposits=[],
	mcAs = None,			      # Replicate MC match as the one used for Muons
	);				      #  you can specify more than one collection for this

removeMCMatching(process, ['All'])
#removeMCMatching(process, ['All'],outputInProcess = False)


# Bs->JPsif0 reconstruction done in module with coded fit
process.mkcands = cms.EDAnalyzer("Bfinder",
  VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
  HLTriggerResults = cms.untracked.string('HLT'),
#  HLTriggerResults = cms.untracked.string('REDIGI36X'),
#  GenParticles = cms.untracked.string('genParticlesPlusSim'),
  GenParticles = cms.untracked.string('genParticles'),
  doMC = cms.untracked.bool(False)
)



###################################################################
###################################################################
# New (easier) Onia2mumu trigger matching
#
#    # Make PAT Muons

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()


from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

#
#################################################################
#################################################################


### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""),
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Bfinder_2012X.root')
)


process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
#process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)


# can I do a replace of patMuons with the sequence that includes the trigger matching?
process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)

process.pat = cms.Path( process.patDefaultSequence )

#print(process.pat)

#patAODTrackCandsUnfiltered*patAODTrackCands*electronMatch*(patTrackCandsMCMatch*patTrackCands+patElectrons)+
#muonMatch*patMuons+pfPileUp+pfNoPileUp*(pfAllNeutralHadrons+pfAllChargedHadrons+pfAllPhotons)*tauIsoDepositPFCandidates*
#tauIsoDepositPFChargedHadrons*tauIsoDepositPFNeutralHadrons*tauIsoDepositPFGammas*tauMatch*tauGenJets*
#tauGenJetsSelectorAllHadrons*tauGenJetMatch*patTaus+photonMatch*patPhotons+patCandidateSummary*selectedPatElectrons+
#selectedPatTrackCands+selectedPatMuons+selectedPatTaus+selectedPatPhotons+selectedPatCandidateSummary*
#cleanPatMuons*(cleanPatElectrons+cleanPatTrackCands)*cleanPatPhotons*cleanPatTaus*cleanPatCandidateSummary*
#countPatElectrons+countPatMuons+countPatTaus+countPatLeptons+countPatPhotons


process.MessageLogger.suppressWarning.append('patTrigger')
process.MessageLogger.cerr.FwkJob.limit=1
process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )


process.ntup = cms.Path(process.mkcands )
process.filter = cms.Path(process.hltLevel1GTSeed * process.noScraping)
process.schedule = cms.Schedule(process.filter, process.pat, process.ntup )


#def customise(process):

process.MessageLogger.cerr.FwkReport.reportEvery = 501

#return (process)
#process = customise(process)
