import FWCore.ParameterSet.Config as cms
from math import pi
from L1Trigger.L1CaloTrigger.Phase1L1TSumsProducer_cfi import Phase1L1TSumsProducer
from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring("file:/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719_oldhgc.done/TTbar_PU200/inputs104X_1.root"),
)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('ComputeJetsAndSums_TTBar_PU200_104XMTD.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer_*_*",
    "keep *_Phase1PFL1TJetProducer_*_*",
    "keep *_ak4GenJetsNoNu__*",
    "keep *_Phase1L1TJetCalibrator_*_*",
    "keep *_Phase1L1TSumsProducer_*_*",
    "keep *_Phase1PFL1TSumsProducer_*_*",
  ),
)

#eta binning of the HLS FW
caloEtaSegmentation = cms.vdouble(
  0.0, 0.0833, 0.1666, 0.2499, 0.3332, 0.4165, 0.4998, 0.5831, 0.6664, 0.7497, 
  0.833, 0.9163, 0.9996, 1.0829, 1.1662, 1.2495, 1.3328, 1.4161, 1.5
)

# sets up jet finder
process.Phase1L1TJetProducer = cms.EDProducer('Phase1L1TJetProducer',
  inputCollectionTag = cms.InputTag("l1pfCandidates", "Puppi"),
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(8),
  phiLow = cms.double(0),
  phiUp = cms.double(0.7),
  jetIEtaSize = cms.uint32(5),
  jetIPhiSize = cms.uint32(5),
  seedPtThreshold = cms.double(5), # GeV
  puSubtraction = cms.bool(False),
  outputCollectionName = cms.string("UncalibratedPhase1L1TJetFromPfCandidates"),
  vetoZeroPt = cms.bool(True)
)

# lut configuration, you can generate your own with test/generateConSinPhi.py 
sinPhi = cms.vdouble(0.04374, 0.13087, 0.21701, 0.30149, 0.38365, 0.46289, 0.53858, 0.61015)
cosPhi = cms.vdouble(0.99904, 0.9914, 0.97617, 0.95347, 0.92348, 0.88642, 0.84257, 0.79229)

# sum module
process.Phase1L1TSumsProducer = cms.EDProducer('Phase1L1TSumsProducer',
  particleCollectionTag = cms.InputTag("l1pfCandidates", "Puppi"),
  jetCollectionTag = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidates"),
  nBinsPhi = cms.uint32(8),
  phiLow = cms.double(0),
  phiUp = cms.double(0.7),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  htPtThreshold = cms.double(30),
  outputCollectionName = cms.string("Sums"),
)

# runs the jet finder and sum producer
process.p = cms.Path(
  process.Phase1L1TJetProducer + 
  process.Phase1L1TSumsProducer
)

process.e = cms.EndPath(process.out)