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

#include <cmath>

class Phase1L1TSumsProducer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit Phase1L1TSumsProducer(const edm::ParameterSet&);
      ~Phase1L1TSumsProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      l1t::EtSum _computeHT(const std::vector<reco::CaloJet>& l1jetVector);
      l1t::EtSum _computeMHT(const std::vector<reco::CaloJet>& l1jetVector);

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
      // lower eta value
      double _etaLow;
      // higher eta value
      double _etaUp;
      // size of a phi bin
      double _phiStep;
      // threshold for ht calculation
      double _htPtThreshold;
      // threshold for ht calculation
      double _mhtPtThreshold;
      // jet eta cut for ht calculation
      double _htAbsEtaCut;
      // jet eta cut for mht calculation
      double _mhtAbsEtaCut;
      // label of sums
      std::string _outputCollectionName;

      bool _debug;

};

// initialises plugin configuration and prepares ROOT file for saving the sums
Phase1L1TSumsProducer::Phase1L1TSumsProducer(const edm::ParameterSet& iConfig):
  // getting configuration settings
  _sinPhi(iConfig.getParameter<std::vector<double> >("sinPhi")),
  _cosPhi(iConfig.getParameter<std::vector<double> >("cosPhi")),
  _nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
  _phiLow(iConfig.getParameter<double>("phiLow")),
  _phiUp(iConfig.getParameter<double>("phiUp")),
  _etaLow(iConfig.getParameter<double>("etaLow")),
  _etaUp(iConfig.getParameter<double>("etaUp")),
  _htPtThreshold(iConfig.getParameter<double>("htPtThreshold")),
  _mhtPtThreshold(iConfig.getParameter<double>("mhtPtThreshold")),
  _htAbsEtaCut(iConfig.getParameter<double>("htAbsEtaCut")),
  _mhtAbsEtaCut(iConfig.getParameter<double>("mhtAbsEtaCut")),
  _outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName")),
  _debug(iConfig.getParameter<bool>("debug"))
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
  produces< std::vector<l1t::EtSum> >( this -> _outputCollectionName ).setBranchAlias(this -> _outputCollectionName);
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
  l1t::EtSum lMHT = this -> _computeMHT(*jetCollectionHandle);

  //packing sums in vector for event saving
  std::unique_ptr< std::vector<l1t::EtSum> > lSumVectorPtr(new std::vector<l1t::EtSum>(0));
  lSumVectorPtr -> push_back(lHT);
  lSumVectorPtr -> push_back(lMHT);
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
  double lHT = 0;
  // range-based for loop that goes through all the trigger jets in the event
  for (const auto & jet: l1jetVector)
  {
    double lJetPt = jet.pt();
    double lJetPhi = jet.phi();
    double lJetEta = jet.eta();
    if 
    (
      (lJetPhi < this -> _phiLow) ||
      (lJetPhi >= this -> _phiUp)  ||
      (lJetEta < this -> _etaLow)||
      (lJetEta >= this -> _etaUp)
    ) continue;

    // std::cout << "HT : " << lJetPt >= this -> _htPtThreshold << " " << std::fabs( lJetEta ) < _htAbsEtaCut << " " << (lJetPt >= this -> _htPtThreshold && std::fabs( lJetEta ) < _htAbsEtaCut ) << " " << lHT << std::endl;
    lHT += (lJetPt >= this -> _htPtThreshold && std::fabs( lJetEta ) < _htAbsEtaCut ) ? lJetPt : 0;
  }

  reco::Candidate::PolarLorentzVector lHTVector;
  lHTVector.SetPt(lHT);
  lHTVector.SetEta(0);
  lHTVector.SetPhi(0);
  // kTotalHt the the EtSum enumerator type for the HT
  // Not setting hwPt etc. yet
  l1t::EtSum lHTSum(lHTVector, l1t::EtSum::EtSumType::kTotalHt, 0, 0, 0, 0);
  return lHTSum;
}

// computes MHT
// adds jet pt to mht only if the pt of the jet is above the mht calculation threshold
l1t::EtSum Phase1L1TSumsProducer::_computeMHT(const std::vector<reco::CaloJet>& l1jetVector)
{
  
  unsigned int lTotalJetPx = 0;
  unsigned int lTotalJetPy = 0;

  std::vector<unsigned int> jetPtInPhiBins(_nBinsPhi, 0);

  for (const auto & jet: l1jetVector)
  { 
    double lJetPhi = jet.phi();

    if ((lJetPhi < this -> _phiLow) || (lJetPhi >= this -> _phiUp)) continue;

    unsigned int iPhi = ( lJetPhi - this -> _phiLow ) / this -> _phiStep;

    if ( jet.pt() >= this -> _mhtPtThreshold && std::fabs( jet.eta() ) < _mhtAbsEtaCut ) {
      unsigned int digiJetPt = floor( jet.pt() / 0.25 );
      jetPtInPhiBins[iPhi] += digiJetPt;

      if ( _debug ) {
        std::cout << "Jet pt : " << jet.pt() << " " << jet.pt() / 0.25 << " " << jet.eta() << " " << jet.phi() << std::endl;
      }
    }
  }

  if ( _debug ) {
    std::cout << "Jet pt sums in phi bins : " << std::endl;
    for (const auto& pt : jetPtInPhiBins ) {
      std::cout << pt << " " << std::endl;
    }
  }

  for ( unsigned int iPhi = 0; iPhi < jetPtInPhiBins.size(); ++iPhi ) {

    unsigned int digiJetPtSum = jetPtInPhiBins[iPhi];

    // retrieving sin cos from LUT emulator
    double lSinPhi = this -> _sinPhi[iPhi];
    double lCosPhi = this -> _cosPhi[iPhi];

    if ( _debug ) {
      if ( digiJetPtSum > 0 ) {
        std::cout << "Jet pt sum : " << digiJetPtSum << std::endl;
        std::cout << "sin cos : " << lSinPhi << " " << lCosPhi << std::endl;
        std::cout << "Px, py : " << floor( digiJetPtSum * lCosPhi ) << " " << floor( digiJetPtSum * lSinPhi ) << std::endl;
        std::cout << "lTotalJetPx/Py : " << lTotalJetPx << " " << lTotalJetPy << std::endl;
      }
    }
    
    // checking if above threshold
    lTotalJetPx += floor( digiJetPtSum * lCosPhi );
    lTotalJetPy += floor( digiJetPtSum * lSinPhi );

  }

  double lMHT = floor( sqrt(lTotalJetPx * lTotalJetPx + lTotalJetPy * lTotalJetPy) ) * 0.25;

  // math::XYZTLorentzVector lMHTVector( lTotalJetPx, lTotalJetPy, 0, digi_lMHT * 0.25 );
  math::PtEtaPhiMLorentzVector lMHTVector( lMHT, 0, acos(lTotalJetPx / ( lMHT / 0.25 )), 0 );
  l1t::EtSum lMHTSum(lMHTVector, l1t::EtSum::EtSumType::kMissingHt, 0, 0, 0, 0 );
  // kTotalMHT the the EtSum enumerator type for the MHT
  // l1t::EtSum lMHTSum(lMHTVector, l1t::EtSum::EtSumType::kMissingHt, 0, 0, 0, 0 );

  if ( _debug ) std::cout << "Output MHT : " << lMHT << " " << lMHTSum.pt() << std::endl;
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
