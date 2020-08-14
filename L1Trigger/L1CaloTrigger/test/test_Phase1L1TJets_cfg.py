import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")
process.TFileService = cms.Service('TFileService', fileName = cms.string("MyTTrees.root"))
process.SaveSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
  genMETCollectionTag = cms.InputTag("genMetTrue"), # taking pre-existing MET collection
  l1tMETCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"), # taking L1T MET produced by jet trigger
  genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"), # taking pre-existing gen jet collection
  l1tHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"), # taking L1T HT produced by jet trigger
  l1tMHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"), # taking L1T MHT produced by jet trigger

)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

# fileList = FileUtils.loadListFromFile('ttbar.list')
# readFiles = cms.untracked.vstring(*fileList)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/QCD_TTBar_PU200_10_4_0_MTD/TTBAR_QCD_Merged.root"),
  #fileNames = cms.untracked.vstring("file:/hdfs/user/sb17498/CMS_Phase_2/jetMETStudies/SingleNeutrino_PU200_104XMTD/SingleNeutrino_PU200.root"),
  #fileNames = cms.untracked.vstring(
  #  "file:pf500.root",
  #)
  # skipEvents = cms.untracked.uint32(181)
)

# Loads 7x7 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

# Load 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9_cff')

# Load trimmed 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9trimmed_cff')


# AK4 PF jets
process.load('L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff')
process.l1PFJets = cms.Sequence( process.ak4PFL1Puppi + process.ak4PFL1PuppiCorrected )

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile_debug.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer*_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetCalibrator*_*_*",
    "keep *_Phase1L1TSumsProducer*_*_*",
    "keep *_ak4PFL1Puppi*_*_*",
    "keep *_l1PFMetPuppi*_*_*",
    "keep *_genMetTrue_*_*"
    # "keep nanoaodFlatTable_*Table_*_*"
  ),
)

# process.SaveSums = cms.EDProducer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("", "", "")
#   l1tMETCollectionTag = cms.InputTag("", "", "")
#   genJetCollectionTag = cms.InputTag("", "", "")
#   l1tHTCollectionTag = cms.InputTag("", "", "")
# )

# process.l1pfjetTable = cms.EDProducer("L1PFJetTableProducer",
#     gen = cms.InputTag("ak4GenJetsNoNu"),
#     commonSel = cms.string("pt > 5 && abs(eta) < 5.0"),
#     drMax = cms.double(0.2),
#     minRecoPtOverGenPt = cms.double(0.1),
#     jets = cms.PSet(
#         Gen = cms.InputTag("ak4GenJetsNoNu"),
#         Gen_sel = cms.string("pt > 15"),
#         AK4 = cms.InputTag("ak4PFL1Puppi"),
#         # SeedCone = cms.InputTag("scL1Puppi"),
#         PhaseI7x7 = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidates"),
#         PhaseI9x9 = cms.InputTag("Phase1L1TJetProducer9x9", "UncalibratedPhase1L1TJetFromPfCandidates"),
#         PhaseI9x9trimmed = cms.InputTag("Phase1L1TJetProducer9x9trimmed", "UncalibratedPhase1L1TJetFromPfCandidates")
#     ),
#     moreVariables = cms.PSet(
#     ),
# )

# process.outnano = cms.OutputModule("NanoAODOutputModule",
#     fileName = cms.untracked.string("perfNano.root"),
#     SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
#     outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
#     compressionLevel = cms.untracked.int32(4),
#     compressionAlgorithm = cms.untracked.string("ZLIB"),
# )

# process.p = cms.Path(process.Phase1L1TJetsSequence * process.Phase1L1TJetsSequence9x9 * process.Phase1L1TJetsSequence9x9trimmed * process.l1PFJets * process.l1PFMetPuppi * process.l1pfjetTable )
process.p = cms.Path(process.Phase1L1TJetsSequence * process.Phase1L1TJetsSequence9x9 * process.Phase1L1TJetsSequence9x9trimmed * process.l1PFJets * process.l1PFMetPuppi )
# process.p = cms.Path(process.Phase1L1TJetsSequence )

process.e = cms.EndPath(process.out * process.outnano)


# process.out = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string('myOutputFile.root'),
#   outputCommands = cms.untracked.vstring(
#     "drop *",
#     "keep *_Phase1L1TJetProducer_*_*",
#     "keep *_ak4GenJetsNoNu_*_*",
#     "keep *_Phase1L1TJetCalibrator_*_*",
#   ),
# )

# from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer

# Phase1L1TJetsSequence = cms.Sequence(
#   Phase1L1TJetProducer
# )

# process.p = cms.Path(process.Phase1L1TJetsSequence)

# process.e = cms.EndPath(process.out)