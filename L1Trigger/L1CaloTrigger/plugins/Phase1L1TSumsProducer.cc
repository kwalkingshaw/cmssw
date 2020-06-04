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
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include <cmath>
#include "ap_fixed.h"
#include "ap_int.h"

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
      l1t::EtSum _computeMHT(const std::vector<reco::CaloJet>& l1jetVector);

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_particleCollectionTag;
      edm::EDGetTokenT<std::vector<reco::CaloJet> > *_jetCollectionTag;

      // holds the sin and cos for HLs LUT emulation
      ap_ufixed<8, 1, AP_RND> _sinPhi; 
      ap_ufixed<8, 1, AP_RND> _cosPhi; 
      unsigned int _nBinsPhi;
      
      double _lsb;
      double _lsb_pt;
      // lower phi value
      ap_uint<8> _phiLow_hls;
      // higher phi value
      ap_uint<8> _phiUp_hls;
      // lower eta value
      ap_uint<8> _etaLow_hls;
      // higher eta value
      ap_uint<8> _etaUp_hls;
      // size of a phi bin
      ap_uint<8> _phiStep_hls;
      // threshold for ht calculation
      ap_uint<16>_htPtThreshold_hls;
      // label of sums
      std::string _outputCollectionName;

    

};

// initialises plugin configuration and prepares ROOT file for saving the sums
Phase1L1TSumsProducer::Phase1L1TSumsProducer(const edm::ParameterSet& iConfig):
  // getting configuration settings
  _lsb(iConfig.getParameter<double>("lsb")),
  _lsb_pt(iConfig.getParameter<double>("lsb_pt")),
  _sinPhi(iConfig.getParameter<std::vector< ap_ufixed<8, 1, AP_RND> >>("sinPhi")),
  _cosPhi(iConfig.getParameter<std::vector< ap_ufixed<8, 1, AP_RND> >>("cosPhi")),
  _nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi") / _lsb),
  _phiLow_hls(iConfig.getParameter<double>("phiLow") / _lsb),
  _phiUp_hls(iConfig.getParameter<double>("phiUp") / _lsb),
  _etaLow_hls(iConfig.getParameter<double>("etaLow") / _lsb),
  _etaUp_hls(iConfig.getParameter<double>("etaUp") / _lsb),
  _htPtThreshold_hls(iConfig.getParameter<double>("htPtThreshold") / _lsb_pt),
  _outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName"))
  

{
  // three things are happening in this line:
  // 1) retrieving the tag for the input particle collection with "iConfig.getParameter(string)"
  // 2) telling CMSSW that I will retrieve a collection of pf candidates later with "consumes< edm::View<reco::Candidate>(InputTag)"
  // 3) obtaining a token that will enable me to access data later with "new edm::EDGetTokenT< edm::View<reco::Candidate> >""
  this -> _particleCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("particleCollectionTag")));  
  // same thing here, I am setting myself up to access jets down the road
  this -> _jetCollectionTag = new edm::EDGetTokenT< std::vector<reco::CaloJet> >(consumes< std::vector<reco::CaloJet> > (iConfig.getParameter< edm::InputTag >("jetCollectionTag")));  
  this -> _phiStep_hls = ( this -> _phiUp_hls - this -> _phiLow_hls ) / this -> _nBinsPhi;
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
  l1t::EtSum lMHT = this -> _computeMHT(*jetCollectionHandle);

  //packing sums in BXVector for event saving
  std::unique_ptr< BXVector<l1t::EtSum> > lSumVectorPtr(new BXVector<l1t::EtSum>(2));
  lSumVectorPtr -> push_back(0, lHT);
  lSumVectorPtr -> push_back(0, lMET);
  lSumVectorPtr -> push_back(0, lMHT);
  //std::cout << "HT-MET sums prod: " << lHT.pt() << "\t" << lMET.pt() << std::endl;
  //std::cout << "MET-MHT sums prod: " << lMET.pt() << "\t" << lMHT.pt() << std::endl;
  //std::cout << "MHT sums prod: " << lMHT.pt() << std::endl;
  //saving sums
  iEvent.put(std::move(lSumVectorPtr), this -> _outputCollectionName);

  return;

}

// computes ht, adds jet pt to ht only if the pt of the jet is above the ht calculation threshold
l1t::EtSum Phase1L1TSumsProducer::_computeHT(const std::vector<reco::CaloJet>& l1jetVector) 
{
  ap_uint<16> lHT = 0;
  
  // range-based for loop that goes through all the trigger jets in the event
  for (const auto & jet: l1jetVector)
  {
    ap_uint<16> lJetPt = jet.pt();
    ap_uint<8> lJetPhi = jet.phi();
    ap_uint<8> lJetEta = jet.eta();
    if 
    (
      (lJetPhi < this -> _phiLow_hls) ||
      (lJetPhi >= this -> _phiUp_hls)  ||
      (lJetEta < this -> _etaLow_hls)||
      (lJetEta >= this -> _etaUp_hls)  
    ) continue;
    //lHT += (lJetPt >= this -> _htPtThreshold_hls) ? lJetPt : 0;
    //lHT *= lsb_pt;
    if (lJetPt >= this -> _htPtThreshold_hls){
      lHT = lHT + (lJetPt * _lsb_pt); 
      //lHT += lJetPt * lsb_pt;
      
    }
    
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

  ap_int<16> lTotalPx = 0;
  ap_int<16> lTotalPy = 0;
  
  // range-based for loop, that goes through the particle flow inputs
  for (const auto & particle : particleCollection)
  {
    ap_uint<8> lParticlePhi = particle.phi();
    ap_uint<8> lParticleEta = particle.eta();
    ap_uint<16> lParticlePt = particle.pt();
    // skip particle if it does not fall in the range
    if 
    (
      (lParticlePhi < this -> _phiLow_hls) ||
      (lParticlePhi >= this -> _phiUp_hls)  ||
      (lParticleEta < this -> _etaLow_hls)  ||
      (lParticleEta >= this -> _etaUp_hls)  
    ) continue;
    // computing bin index
    unsigned int iPhi = ( lParticlePhi - this -> _phiLow_hls ) / this -> _phiStep_hls;
    // retrieving sin cos from LUT emulator
    ap_ufixed<8, 1, AP_RND> lSinPhi = this -> _sinPhi[iPhi];
    ap_ufixed<8, 1, AP_RND> lCosPhi = this -> _cosPhi[iPhi];
    // computing px and py of the particle and adding it to the total px and py of the event
    lTotalPx += (lParticlePt * lCosPhi);
    lTotalPy += (lParticlePt * lSinPhi);
  }

  ap_uint<16> lMET = sqrt(static_cast<double>(lTotalPx * lTotalPx + lTotalPy * lTotalPy)) * _lsb_pt;
  //packing in EtSum object
  reco::Candidate::PolarLorentzVector lMETVector;
  double lCosMETPhi = lTotalPx/lMET;
  lMETVector.SetPxPyPzE(lTotalPx, lTotalPy, 0, lMET);
  // lMETVector.SetEta(0);
  // lMETVector.SetPhi(0);
  // kMissingEt is the enumerator for MET
  l1t::EtSum lMETSum(lMETVector, l1t::EtSum::EtSumType::kMissingEt);
  // std::cout << lMETVector.pt() << " == " << lMET << "?" << std::endl;

  return lMETSum;

}

// computes MHT
// adds jet pt to mht only if the pt of the jet is above the mht calculation threshold
l1t::EtSum Phase1L1TSumsProducer::_computeMHT(const std::vector<reco::CaloJet>& l1jetVector)
{
  
  ap_int<16> lTotalJetPx = 0;
  ap_int<16> lTotalJetPy = 0;
  
  

  for (const auto & jet: l1jetVector)
  { 
    ap_uint<16> lJetPhi = jet.phi();
    ap_uint<16> lJetPt = jet.pt();
    

    if ((lJetPhi < this -> _phiLow_hls) || (lJetPhi >= this -> _phiUp_hls)) continue;

    unsigned int iPhi = ( lJetPhi - this -> _phiLow_hls ) / this -> _phiStep_hls;

    // retrieving sin cos from LUT emulator
    ap_ufixed<8, 1, AP_RND> lSinPhi = this -> _sinPhi[iPhi];
    ap_ufixed<8, 1, AP_RND> lCosPhi = this -> _cosPhi[iPhi];
    
   
    // checking if above threshold
    //lTotalJetPx += (jet.pt() >= this -> _htPtThreshold_hls) ? jet.pt() * lCosPhi : 0;
    //lTotalJetPy += (jet.pt() >= this -> _htPtThreshold_hls) ? jet.pt() * lSinPhi : 0;

    if (lJetPt >= this -> _htPtThreshold_hls){
      lTotalJetPx += lJetPt * lCosPhi;
      lTotalJetPy += lJetPt * lSinPhi;
      
    }

  
  }

  ap_uint<16> lMHT = sqrt(static_cast<double> (lTotalJetPx * lTotalJetPx + lTotalJetPy * lTotalJetPy)) * _lsb_pt;

  reco::Candidate::PolarLorentzVector lMHTVector;
  lMHTVector.SetPt(lMHT);
  lMHTVector.SetEta(0);
  lMHTVector.SetPhi(0);

  // kTotalMHT the the EtSum enumerator type for the MHT
  l1t::EtSum lMHTSum(lMHTVector, l1t::EtSum::EtSumType::kMissingHt);

  return lMHTSum;

}

// I have no idea what this does, it is created by default by CMSSW
void Phase1L1TSumsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// creates the plugin for later use in python
DEFINE_FWK_MODULE(Phase1L1TSumsProducer);
