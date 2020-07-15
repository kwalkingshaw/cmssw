import FWCore.ParameterSet.Config as cms
from math import pi

sinPhi = cms.vdouble(0.04374, 0.13087, 0.21701, 0.30149, 0.38365, 0.46289, 0.53858, 0.61015, 0.67705, 0.73877,
                     0.79484, 0.84483, 0.88835, 0.92508, 0.95473, 0.97707, 0.99194, 0.99922, 0.99885, 0.99084,
                     0.97525, 0.9522, 0.92186, 0.88446, 0.8403, 0.78971, 0.73308, 0.67084, 0.60347, 0.53148,
                     0.45542, 0.37588, 0.29346, 0.2088, 0.12253, 0.03534, -0.05213, -0.1392, -0.22521, -0.30949,
                     -0.3914, -0.47032, -0.54564, -0.61679, -0.68322, -0.74441, -0.79992, -0.8493, -0.89218,
                     -0.92824, -0.9572, -0.97883, -0.99297, -0.99952, -0.99841, -0.98967, -0.97336, -0.94959,
                     -0.91857, -0.88051, -0.83572, -0.78453, -0.72734, -0.66458, -0.59674, -0.52434, -0.44792,
                     -0.36807, -0.28541, -0.20057, -0.11419, -0.02693)
cosPhi = cms.vdouble(0.99904, 0.9914, 0.97617, 0.95347, 0.92348, 0.88642, 0.84257, 0.79229, 0.73593, 0.67395,
                     0.60681, 0.53503, 0.45916, 0.37977, 0.29747, 0.2129, 0.1267, 0.03954, -0.04794, -0.13504,
                     -0.22111, -0.30549, -0.38753, -0.46661, -0.54212, -0.61348, -0.68014, -0.7416, -0.79739, -0.84707,
                     -0.89028, -0.92667, -0.95597, -0.97796, -0.99246, -0.99938, -0.99864, -0.99026, -0.97431, -0.9509,
                     -0.92022, -0.88249, -0.83802, -0.78713, -0.73022, -0.66772, -0.60011, -0.52791, -0.45167,
                     -0.37198, -0.28944, -0.20468, -0.11836, -0.03113, 0.05633, 0.14337, 0.2293, 0.31349,
                     0.39527, 0.47403, 0.54916, 0.62009, 0.68628, 0.74721, 0.80243, 0.85151, 0.89407,
                     0.9298, 0.95841, 0.97968, 0.99346, 0.99964)

#caloEtaSegmentation = cms.vdouble(
#    0.0, 0.0833, 0.1666, 0.2499, 0.3332, 0.4165, 0.4998, 0.5831, 0.6664, 0.7497, 
#    0.833, 0.9163, 0.9996, 1.0829, 1.1662, 1.2495, 1.3328, 1.4161, 1.5)
# sinPhi = cms.vdouble(0.04374, 0.13087, 0.21701, 0.30149, 0.38365, 0.46289, 0.53858, 0.61015),
# cosPhi = cms.vdouble(0.99904, 0.9914, 0.97617, 0.95347, 0.92348, 0.88642, 0.84257, 0.79229),


Phase1L1TSumsProducer = cms.EDProducer('Phase1L1TSumsProducer',
  particleCollectionTag = cms.InputTag("l1pfCandidates", "Puppi"),
  jetCollectionTag = cms.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates"),
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.15),
  phiUp = cms.double(3.15),
  # nBinsPhi = cms.uint32(8),
  # phiLow = cms.double(0),
  # phiUp = cms.double(0.7),
  etaLow = cms.double(-5),
  etaUp = cms.double(5),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  htPtThreshold = cms.double(30),
  htAbsEtaCut = cms.double(2.4),
  mhtPtThreshold = cms.double(30),
  mhtAbsEtaCut = cms.double(2.4),
  outputCollectionName = cms.string("Sums"),
  debug = cms.bool(False)
)
