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
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer*_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetCalibrator*_*_*",
    "keep *_Phase1L1TSumsProducer*_*_*",
    "keep *_ak4PFL1Puppi*_*_*",
    "keep *_l1PFMetPuppi*_*_*",
    "keep *_genMetTrue_*_*",
  ),
)

# process.SaveSums = cms.EDProducer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("", "", "")
#   l1tMETCollectionTag = cms.InputTag("", "", "")
#   genJetCollectionTag = cms.InputTag("", "", "")
#   l1tHTCollectionTag = cms.InputTag("", "", "")
# )

process.p = cms.Path(process.Phase1L1TJetsSequence * process.Phase1L1TJetsSequence9x9 * process.Phase1L1TJetsSequence9x9trimmed * process.l1PFJets * process.l1PFMetPuppi )

process.e = cms.EndPath(process.out)


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