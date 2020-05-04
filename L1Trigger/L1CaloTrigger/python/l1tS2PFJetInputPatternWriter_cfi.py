
import FWCore.ParameterSet.Config as cms

l1tS2PFJetInputPatternWriter = cms.EDAnalyzer('L1TS2PFJetInputPatternWriter',
    pfTag = cms.InputTag("l1pfCandidates", "Puppi", "IN"),
    filename       = cms.untracked.string("pattern"),
    outDir       = cms.untracked.string("patterns"),
    # nChanPerQuad   = cms.untracked.uint32(4),
    # nQuads         = cms.untracked.uint32(18),
    nPayloadFrames = cms.untracked.uint32(15),
    framesPerFile = cms.untracked.uint32(1015),
    evPerFile = cms.untracked.uint32(50)
)


# 50 * 15 = 750, leaving ~ 250 frames to finish processing
