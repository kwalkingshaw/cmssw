
import FWCore.ParameterSet.Config as cms

l1tS2PFJetInputPatternWriter = cms.EDAnalyzer('L1TS2PFJetInputPatternWriter',
    pfTag = cms.InputTag("l1pfCandidates", "Puppi", "IN"),
    filename       = cms.untracked.string("pattern.txt"),
    outDir       = cms.untracked.string("patterns"),
    # nChanPerQuad   = cms.untracked.uint32(4),
    # nQuads         = cms.untracked.uint32(18),
    nHeaderFrames  = cms.untracked.uint32(0),
    nPayloadFrames = cms.untracked.uint32(15),
    nClearFrames   = cms.untracked.uint32(0)
)
