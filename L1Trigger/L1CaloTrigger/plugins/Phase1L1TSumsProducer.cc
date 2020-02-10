// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TSumsProducer
//
/**\class Phase1L1TSumsProducer Phase1L1TSumsProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TSumsProducer.cc

Description: Produces jets with a phase-1 like sliding window algorithm using a collection of reco::Candidates in input.
  If flag is enabled, computes MET and MHT.

*** INPUT PARAMETERS ***
  * etaBinning: vdouble with eta binning (allows non-homogeneous binning in eta)
  * nBinsPhi: uint32, number of bins in phi
  * phiLow: double, min phi (typically -pi)
  * phiUp: double, max phi (typically +pi)
  * jetIEtaSize: uint32, jet cluster size in ieta
  * jetIPhiSize: uint32, jet cluster size in iphi
  * seedPtThreshold: double, threshold of the seed tower
  * puSubtraction: bool, runs chunky doughnut pile-up subtraction, 9x9 jet only
  * outputCollectionName: string, tag for the output collection
  * vetoZeroPt: bool, controls whether jets with 0 pt should be save. 
    It matters if PU is ON, as you can get negative or zero pt jets after it.
  * outputSumsCollectionName: string, tag for the output collection containing MET and HT
  * inputCollectionTag: inputtag, collection of reco::candidates used as input to the algo
  * htPtThreshold: double, threshold for computing ht

*/
//
// Original Simone Bologna
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

class Phase1L1TSumsProducer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit Phase1L1TSumsProducer(const edm::ParameterSet&);
      ~Phase1L1TSumsProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      l1t::EtSum _computeHT(const std::vector<reco::CaloJet>& l1jetVector);
      template <class ParticleCollection>
      l1t::EtSum _computeMET(const ParticleCollection & caloGrid);

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_particleCollectionTag;
      edm::EDGetTokenT<std::vector<reco::CaloJet> > *_jetCollectionTag;

      // holds the sin and cos for HLs LUT emulation
      std::vector<double> _sinPhi;
      std::vector<double> _cosPhi;
      unsigned int _nBinsPhi;
      // lower phi value
      double _phiLow;
      // higher phi value
      double _phiUp;
      // size of a phi bin
      double _phiStep;
      // threshold for ht calculation
      double _htPtThreshold;
      // label of sums
      std::string _outputCollectionName;

};

// initialises plugin configuration and prepares ROOT file for saving the sums
Phase1L1TSumsProducer::Phase1L1TSumsProducer(const edm::ParameterSet& iConfig):
  // getting configuration settings
  _sinPhi(iConfig.getParameter<std::vector<double> >("sinPhi")),
  _cosPhi(iConfig.getParameter<std::vector<double> >("cosPhi")),
  _nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
  _phiLow(iConfig.getParameter<double>("phiLow")),
  _phiUp(iConfig.getParameter<double>("phiUp")),
  _htPtThreshold(iConfig.getParameter<double>("htPtThreshold")),
  _outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName"))
{
  // three things are happening in this line:
  // 1) retrieving the tag for the input particle collection with "iConfig.getParameter(string)"
  // 2) telling CMSSW that I will retrieve a collection of pf candidates later with "consumes< edm::View<reco::Candidate>(InputTag)"
  // 3) obtaining a token that will enable me to access data later with "new edm::EDGetTokenT< edm::View<reco::Candidate> >""
  this -> _particleCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("particleCollectionTag")));  
  // same thing here, I am setting myself up to access jets down the road
  this -> _jetCollectionTag = new edm::EDGetTokenT< std::vector<reco::CaloJet> >(consumes< std::vector<reco::CaloJet> > (iConfig.getParameter< edm::InputTag >("jetCollectionTag")));  
  this -> _phiStep = ( this -> _phiUp - this -> _phiLow ) / this -> _nBinsPhi;
  // preparing CMSSW to save my sums later
  // "setBranchAlias" specifies the label that my output will have in the output file
  // produces <> sets up the producer to save stuff later
  produces< BXVector<l1t::EtSum> >( this -> _outputCollectionName ).setBranchAlias(this -> _outputCollectionName);
}

// delete dynamically allocated tags (which were created with new)
Phase1L1TSumsProducer::~Phase1L1TSumsProducer()
{
  delete this -> _particleCollectionTag;
  delete this -> _jetCollectionTag;
}

// creates object to access the jets and the pf candidates collection
// then computes the sums and stores them in the EDM root file
void Phase1L1TSumsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // accessing my data
  // handle is an intermediate object between the data on file and on my local memory
  edm::Handle < edm::View< reco::Candidate > > particleCollectionHandle;
  edm::Handle < std::vector<reco::CaloJet> > jetCollectionHandle;
  // this retrieves my data and puts it up for use in my local memory, the handle now behaves like a pointer to the memory area holding my data
  // i.e. I can access it with the * operator or -> operator
  iEvent.getByToken(*(this -> _particleCollectionTag), particleCollectionHandle);
  iEvent.getByToken(*(this -> _jetCollectionTag), jetCollectionHandle);
  
  // computing sums and storing them in sum object
  l1t::EtSum lHT = this -> _computeHT(*jetCollectionHandle);
  l1t::EtSum lMET = this -> _computeMET<>(*particleCollectionHandle);

  //packing sums in BXVector for event saving
  std::unique_ptr< BXVector<l1t::EtSum> > lSumVectorPtr(new BXVector<l1t::EtSum>(2));
  lSumVectorPtr -> push_back(0, lHT);
  lSumVectorPtr -> push_back(0, lMET);
  // std::cout << "HT-MET sums prod: " << lHT.pt() << "\t" << lMET.pt() << std::endl;
  //saving sums
  iEvent.put(std::move(lSumVectorPtr), this -> _outputCollectionName);

  return;

}

// computes ht, adds jet pt to ht only if the pt of the jet is above the ht calculation threshold
l1t::EtSum Phase1L1TSumsProducer::_computeHT(const std::vector<reco::CaloJet>& l1jetVector) 
{
  double lHT = 0;
  // range-based for loop that goes through all the trigger jets in the event
  for (const auto & jet: l1jetVector)
  {
    lHT += (jet.pt() >= this -> _htPtThreshold) ? jet.pt() : 0;
  }

  reco::Candidate::PolarLorentzVector lHTVector;
  lHTVector.SetPt(lHT);
  lHTVector.SetEta(0);
  lHTVector.SetPhi(0);
  // kTotalHt the the EtSum enumerator type for the HT
  l1t::EtSum lHTSum(lHTVector, l1t::EtSum::EtSumType::kTotalHt);

  return lHTSum;
}

// computes MET
// first computes in which bin of the histogram the input falls in
// the phi bin index is used to retrieve the sin-cos value from the LUT emulator
// the pt of the input is multiplied by that sin cos value to obtain px and py that is added to the total event px & py
// after all the inputs have been processed we compute the total pt of the event, and set that as MET
template <class ParticleCollection>
l1t::EtSum Phase1L1TSumsProducer::_computeMET(const ParticleCollection & particleCollection) 
{

  double lTotalPx = 0;
  double lTotalPy = 0;
  // range-based for loop, that goes through the particle flow inputs
  for (const auto & particle : particleCollection)
  {
    double lParticlePhi = particle.phi();
    double lParticlePt = particle.pt();
    // skip particle if it does not fall in the histogram range
    if ((lParticlePhi < this -> _phiLow) || (lParticlePhi > this -> _phiUp)) continue;
    // computing bin index
    unsigned int iPhi = ( lParticlePhi - this -> _phiLow ) / this -> _phiStep;
    // retrieving sin cos from LUT emulator
    double lSinPhi = this -> _sinPhi[iPhi];
    double lCosPhi = this -> _cosPhi[iPhi];
    // computing px and py of the particle and adding it to the total px and py of the event
    lTotalPx += (lParticlePt * lCosPhi);
    lTotalPy += (lParticlePt * lSinPhi);
  }

  double lMET = sqrt(lTotalPx * lTotalPx + lTotalPy * lTotalPy);
  //packing in EtSum object
  reco::Candidate::PolarLorentzVector lMETVector;
  lMETVector.SetPt(lMET);
  lMETVector.SetEta(0);
  lMETVector.SetPhi(0);
  // kMissingEt is the enumerator for MET
  l1t::EtSum lMETSum(lMETVector, l1t::EtSum::EtSumType::kMissingEt);

  return lMETSum;

}

// I have no idea what this does, it is created by default by CMSSW
void Phase1L1TSumsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// creates the plugin for later use in python
DEFINE_FWK_MODULE(Phase1L1TSumsProducer);
