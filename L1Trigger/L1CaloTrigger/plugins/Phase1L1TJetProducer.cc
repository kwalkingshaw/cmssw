// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TJetProducer
//
/**\class Phase1L1TJetProducer Phase1L1TJetProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TJetProducer.cc

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
// Created: Mon Jul 02 2018
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
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "TH2F.h"

#include <algorithm>

class Phase1L1TJetProducer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit Phase1L1TJetProducer(const edm::ParameterSet&);
      ~Phase1L1TJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      /// Finds the seeds in the caloGrid, seeds are saved in a vector that contain the index in the TH2F of each seed
      std::vector<std::tuple<int, int>> _findSeeds(const TH2F & caloGrid, float seedThreshold);
      
      std::vector<reco::CaloJet> _buildJetsFromSeedsWithPUSubtraction(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds, bool killZeroPt);
      std::vector<reco::CaloJet> _buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds);
      
      void _subtract9x9Pileup(const TH2F & caloGrid, reco::CaloJet & jet);
      
      /// Get the energy of a certain tower while correctly handling phi periodicity in case of overflow
      float _getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi);
      
      reco::CaloJet _buildJetFromSeed(const TH2F & caloGrid, const std::tuple<int, int> & seed);
      
      // <3 handy method to fill the calogrid with whatever type
      template <class Container>
      void _fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives, const unsigned int pfRegionIndex );

      // Arranges the input candidates into containers containing the inputs for each region
      // Might also sort the candidates within each region (e.g. by pt) and truncate the inputs.
      template <class Handle >
      std::vector< std::vector< reco::CandidatePtr > > _prepareInputsIntoRegions( const Handle triggerPrimitives );

      unsigned int pfRegionIndex( const unsigned int phiRegion, const unsigned int etaRegion );
      std::pair< double, double> pfRegionEtaPhiLowEdges( const unsigned int pfRegionIndex );
      std::pair< double, double> pfRegionEtaPhiUpEdges( const unsigned int pfRegionIndex );

      l1t::EtSum _computeMET( const TH2F & caloGrid, double etaCut, l1t::EtSum::EtSumType sumType );

      // Determine if this tower should be trimmed or not
      // Trim == removing 3 towers in each corner of the square grid
      // giving a cross shaped grid, which is a bit more circular in shape than a square
      bool _trimTower( const int etaIndex, const int phiIndex );

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_inputCollectionTag;

      std::vector<double> _etaBinning;
      size_t _nBinsEta;
      unsigned int _nBinsPhi;
      double _phiLow;
      double _phiUp;
      unsigned int _jetIEtaSize;
      unsigned int _jetIPhiSize;
      bool _trimmedGrid;
      double _seedPtThreshold;
      double _philsb;
      double _etalsb;
      bool _puSubtraction;
      bool _vetoZeroPt;
      std::vector< double > _etaPFRegionEdges;
      std::vector< double > _phiPFRegionEdges;
      unsigned int _maxInputsPerPFRegion;
      std::vector<double> _sinPhi;
      std::vector<double> _cosPhi;
      // input eta cut for met calculation
      double _metAbsEtaCut;
      // input eta cut for metHF calculation
      double _metHFAbsEtaCut;
      bool _debug;

      std::string _outputCollectionName;

      // std::vector< std::pair<double, double> > _etaPFRegionEdges;


};


Phase1L1TJetProducer::Phase1L1TJetProducer(const edm::ParameterSet& iConfig):
  // getting configuration settings
  _etaBinning(iConfig.getParameter<std::vector<double> >("etaBinning")),
  _nBinsEta(this -> _etaBinning.size() - 1),
  _nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
  _phiLow(iConfig.getParameter<double>("phiLow")),
  _phiUp(iConfig.getParameter<double>("phiUp")),
  _jetIEtaSize(iConfig.getParameter<unsigned int>("jetIEtaSize")),
  _jetIPhiSize(iConfig.getParameter<unsigned int>("jetIPhiSize")),
  _trimmedGrid(iConfig.getParameter<bool>("trimmedGrid")),
  _seedPtThreshold(iConfig.getParameter<double>("seedPtThreshold")),
  _philsb(iConfig.getParameter<double>("philsb")),
  _etalsb(iConfig.getParameter<double>("etalsb")),
  _puSubtraction(iConfig.getParameter<bool>("puSubtraction")),
  _vetoZeroPt(iConfig.getParameter<bool>("vetoZeroPt")),
  _etaPFRegionEdges(iConfig.getParameter<std::vector<double> >("pfEtaRegions")),
  _phiPFRegionEdges(iConfig.getParameter<std::vector<double> >("pfPhiRegions")),
  _maxInputsPerPFRegion(iConfig.getParameter<unsigned int>("maxInputsPerPFRegion")),
  _sinPhi(iConfig.getParameter<std::vector<double> >("sinPhi")),
  _cosPhi(iConfig.getParameter<std::vector<double> >("cosPhi")),
  _metAbsEtaCut(iConfig.getParameter<double>("metAbsEtaCut")),
  _metHFAbsEtaCut(iConfig.getParameter<double>("metHFAbsEtaCut")),
  _debug(iConfig.getParameter<bool>("debug")),
  _outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName"))
{
  this -> _inputCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("inputCollectionTag")));  
  produces<std::vector<reco::CaloJet> >( this -> _outputCollectionName ).setBranchAlias(this -> _outputCollectionName);
  produces< std::vector<l1t::EtSum> >( this -> _outputCollectionName + "MET" ).setBranchAlias(this -> _outputCollectionName + "MET");

  // Make into an exception...
  if ( !std::is_sorted(_etaPFRegionEdges.begin(),_etaPFRegionEdges.end()) ||
       !std::is_sorted(_phiPFRegionEdges.begin(),_phiPFRegionEdges.end())
       ) {
    std::cout << "Input PF regions are not sorted!!!!" << std::endl;
  }

}

Phase1L1TJetProducer::~Phase1L1TJetProducer()
{
  delete this -> _inputCollectionTag;
}

float Phase1L1TJetProducer::_getTowerEnergy(const TH2F & caloGrid, int iEta, int iPhi)
{
  // We return the pt of a certain bin in the calo grid, taking account of the phi periodicity when overflowing (e.g. phi > phiSize), and returning 0 for the eta out of bounds

  int nBinsEta = caloGrid.GetNbinsX();
  int nBinsPhi = caloGrid.GetNbinsY();
  while (iPhi < 1) {
    iPhi += nBinsPhi;
  }
  while (iPhi > nBinsPhi) {
    iPhi -= nBinsPhi;
  }
  if (iEta < 1) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
  if (iEta > nBinsEta) {
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << 0 << ", " << iEta << ", " << iPhi << ")" << std::endl;
    return 0;
  }
    //std::cout << "Returning (pt, ieta, iphi): " << "(" << caloGrid.GetBinContent(iEta, iPhi) << ", " << iEta << ", " << iPhi << ")" << std::endl;
  return caloGrid.GetBinContent(iEta, iPhi);
}

void Phase1L1TJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  TH2F lCaloGrid("caloGrid", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
  lCaloGrid.GetXaxis() -> SetTitle("#eta");
  lCaloGrid.GetYaxis() -> SetTitle("#phi");

  if ( _debug ) {
    TAxis* phiAxis = lCaloGrid.GetYaxis();
    std::cout << "Phi bins : ";
    for ( int i = 0; i < phiAxis->GetNbins(); ++i ) {
      std::cout << phiAxis->GetBinUpEdge(i) << " ";
    }
    std::cout << std::endl;
    TAxis* etaAxis = lCaloGrid.GetXaxis();
    std::cout << "Eta bins : ";
    for ( int i = 0; i < etaAxis->GetNbins(); ++i ) {
      std::cout << etaAxis->GetBinUpEdge(i) << " ";
    }
    std::cout << std::endl;
  }
  
  if (this -> _inputCollectionTag) {
    edm::Handle < edm::View< reco::Candidate > > inputCollectionHandle;
    iEvent.getByToken(*(this -> _inputCollectionTag), inputCollectionHandle);


    // std::vector< std::pair<double, double> > etaRegionEdges{ {0, 0.75}, {0.75, 1.5} };
    std::vector< std::vector< reco::CandidatePtr > > inputsInRegions = _prepareInputsIntoRegions<>( inputCollectionHandle );

    for (unsigned int i = 0; i < inputsInRegions.size(); ++i ) {
      // for ( const auto& p : inputsInRegions[i] ) {
      //   std::cout << p->pt() << " " << p->eta() << " " << p->phi() << std::endl;
      // }
      this -> _fillCaloGrid<>(lCaloGrid, inputsInRegions[i], i);
    }

    // dumping the data
    // this -> _fillCaloGrid<>(lCaloGrid, *inputCollectionHandle);

    // int nBinsX = lCaloGrid.GetNbinsX();
    // int nBinsY = lCaloGrid.GetNbinsY();
    // for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
    // {
    //   std::cout << "iPhi " << iPhi - 1 << " " << lCaloGrid.GetYaxis() -> GetBinCenter(iPhi) << ": ";
    //   for (int iEta = 1; iEta <= nBinsX; iEta++)
    //   {
    //     std::cout <<lCaloGrid.GetBinContent(iEta, iPhi) << " ";
    //   }
    //   std::cout << std::endl;
    // }
    const auto seedsVector = this -> _findSeeds(lCaloGrid, this -> _seedPtThreshold); // seedPtThreshold = 6
    if ( _debug ) std::cout << "Number of seeds : " << seedsVector.size() << std::endl;
    std::vector<reco::CaloJet> l1jetVector;
    if (this -> _puSubtraction)
    {
      l1jetVector = this -> _buildJetsFromSeedsWithPUSubtraction(lCaloGrid, seedsVector, this -> _vetoZeroPt);
    } else
    {
      l1jetVector = this -> _buildJetsFromSeeds(lCaloGrid, seedsVector);
    }

    //saving jets
    std::unique_ptr< std::vector<reco::CaloJet> > l1jetVectorPtr(new std::vector<reco::CaloJet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), this -> _outputCollectionName);
  }

  if ( _debug ) std::cout << "Calculate MET" << std::endl;

  l1t::EtSum lMET = _computeMET( lCaloGrid, _metAbsEtaCut, l1t::EtSum::EtSumType::kMissingEt );
  l1t::EtSum lMETHF = _computeMET( lCaloGrid, _metHFAbsEtaCut, l1t::EtSum::EtSumType::kMissingEtHF);

  std::unique_ptr< std::vector<l1t::EtSum> > lSumVectorPtr(new std::vector<l1t::EtSum>(0));
  lSumVectorPtr -> push_back(lMET);
  lSumVectorPtr -> push_back(lMETHF);
  // std::cout << "MET sums prod: " << lMET.pt() << std::endl;

  //saving sums
  iEvent.put(std::move(lSumVectorPtr), this -> _outputCollectionName+"MET");


  return;

}

void Phase1L1TJetProducer::_subtract9x9Pileup(const TH2F & caloGrid, reco::CaloJet & jet) {
  //For each jet we compute the 4 side bands
  //for (reco::CaloJet & jet: jetCollection) {
  //for (auto jetIterator = jetCollection.begin(0); jetIterator != jetCollection.end(0); jetIterator++)
  // these variables host the total pt in each sideband and the total pileup contribution
  float topBandPt = 0;
  float leftBandPt = 0;
  float rightBandPt = 0;
  float bottomBandPt = 0;
  float pileUpEnergy;

  // hold the jet's x-y (and z, as I have to use it, even if 2D) location in the histo
  int xCenter, yCenter, zCenter;
  // Retrieving histo-coords for seed
  caloGrid.GetBinXYZ(caloGrid.FindFixBin(jet.eta(), jet.phi()), xCenter, yCenter, zCenter);
  
  // Computing pileup
  for (int x = -4; x <= 4; x++) {
    for (int y = 0; y < 3; y++) {
      // top band, I go up 5 squares to reach the bottom of the top band
      // +x scrolls from left to right, +y scrolls up
      topBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter + (5 + y));
      // left band, I go left 5 squares (-5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls left
      leftBandPt += this -> _getTowerEnergy(caloGrid, xCenter - (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      rightBandPt += this -> _getTowerEnergy(caloGrid, xCenter + (5 + y), yCenter + x);
      // right band, I go right 5 squares (+5) to reach the bottom of the top band
      // +x scrolls from bottom to top, +y scrolls right
      bottomBandPt += this -> _getTowerEnergy(caloGrid, xCenter + x, yCenter - (5 + y));
    }
  }
  // adding bands and removing the maximum band (equivalent to adding the three minimum bands)
  pileUpEnergy = topBandPt + leftBandPt + rightBandPt + bottomBandPt - std::max(topBandPt, std::max(leftBandPt, std::max(rightBandPt, bottomBandPt)));

  //preparing the new 4-momentum vector
  //math::PtEtaPhiMLorentzVector ptVector;
  reco::Candidate::PolarLorentzVector ptVector;
  // removing pu contribution
  float ptAfterPUSubtraction = jet.pt() - pileUpEnergy;
  ptVector.SetPt((ptAfterPUSubtraction > 0)? ptAfterPUSubtraction : 0); 
  ptVector.SetEta(jet.eta());
  ptVector.SetPhi(jet.phi());
  //updating the jet
  jet.setP4(ptVector);
  jet.setPileup(pileUpEnergy);
  return;
}

std::vector<std::tuple<int, int>> Phase1L1TJetProducer::_findSeeds(const TH2F & caloGrid, float seedThreshold) 
{
  int nBinsX = caloGrid.GetNbinsX();
  int nBinsY = caloGrid.GetNbinsY();

  std::vector<std::tuple<int, int>> seeds;

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  // for each point of the grid check if it is a local maximum
  // to do so I take a point, and look if is greater than the points around it (in the 9x9 neighborhood)
  // to prevent mutual exclusion, I check greater or equal for points above and right to the one I am considering (including the top-left point)
  // to prevent mutual exclusion, I check greater for points below and left to the one I am considering (including the bottom-right point)

  for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
  {
    for (int iEta = 1; iEta <= nBinsX; iEta++)
    {
      float centralPt = caloGrid.GetBinContent(iEta, iPhi);
      if (centralPt < seedThreshold) continue;
      bool isLocalMaximum = true;
      if ( _debug ) std::cout << "Testing seed : " << iEta << " " << iPhi << " " << caloGrid.GetXaxis()->GetBinCenter( iEta ) << " " << caloGrid.GetYaxis()->GetBinCenter( iPhi ) << std::endl;
      // Scanning through the grid centered on the seed
      for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
      {
        for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
        {
          if ( _trimmedGrid ) {
            if ( _trimTower( etaIndex, phiIndex ) ) continue;
          }

          if ((etaIndex == 0) && (phiIndex == 0)) continue;
          if (phiIndex > 0) {
            if (phiIndex > -etaIndex){
              isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            } else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            }
          } else {
            if (phiIndex >= -etaIndex){
              isLocalMaximum = ((isLocalMaximum) && (centralPt > this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            } else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex)));
            }
          }
        }
      }
      if (isLocalMaximum)
      {
        if ( _debug ) std::cout << "Seed ieta, iphi : " << iEta << " " << iPhi << std::endl;
        seeds.emplace_back(std::make_tuple(iEta, iPhi));
      }

    }
  }

  return seeds;
}

reco::CaloJet Phase1L1TJetProducer::_buildJetFromSeed(const TH2F & caloGrid, const std::tuple<int, int> & seed) 
{
  int iEta = std::get<0>(seed);
  int iPhi = std::get<1>(seed);

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  if ( _debug ) std::cout << "Building jet based on seed : " << iEta << " " << iPhi << std::endl;

  float ptSum = 0;
  // Scanning through the grid centered on the seed
  for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++)
  {
    for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++)
    {
      if ( _trimmedGrid ) {
        if ( _trimTower( etaIndex, phiIndex ) ) continue;
      }
      if ( _debug ) {
        if ( this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex) > 0 ) {
          std::cout << "Adding this tower to jet : " << this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex) << std::endl;
          std::cout << etaIndex << " " << phiIndex << std::endl;
          std::cout << iEta + etaIndex << " " << iPhi + phiIndex << std::endl;          
        }
      }
      ptSum += this -> _getTowerEnergy(caloGrid, iEta + etaIndex, iPhi + phiIndex);
    }

  }

  // Creating a jet with eta phi centered on the seed and momentum equal to the sum of the pt of the components    
  //reco::Candidate::LorentzVector ptVector;
  reco::Candidate::PolarLorentzVector ptVector;
  ptVector.SetPt(ptSum);
  //ptVector.SetPtEtaPhiE(ptSum, caloGrid.GetXaxis() -> GetBinCenter(iEta), caloGrid.GetYaxis() -> GetBinCenter(iPhi), ptSum);
  ptVector.SetEta(caloGrid.GetXaxis() -> GetBinCenter(iEta));
  ptVector.SetPhi(caloGrid.GetYaxis() -> GetBinCenter(iPhi));
  // ptVector.SetEta(iEta);
  // ptVector.SetPhi(iPhi);
  reco::CaloJet jet;
  jet.setP4(ptVector);
  return jet;
}

std::vector<reco::CaloJet> Phase1L1TJetProducer::_buildJetsFromSeedsWithPUSubtraction(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds, bool killZeroPt)
{

  // For each seed take a grid centered on the seed of the size specified by the user
  // Sum the pf in the grid, that will be the pt of the l1t jet. Eta and phi of the jet is taken from the seed.
  std::vector<reco::CaloJet> jets;
  for (const auto& seed: seeds)
  {
    reco::CaloJet jet = this -> _buildJetFromSeed(caloGrid, seed);
    this -> _subtract9x9Pileup(caloGrid, jet);
    //killing jets with 0 pt
    if ((this -> _vetoZeroPt) && (jet.pt() <= 0)) continue;
    jets.push_back(jet);
  }

  // sort by pt
  std::sort(jets.begin(), jets.end(), 
    [](const reco::CaloJet& jet1, const reco::CaloJet& jet2) 
    {
      return jet1.pt() > jet2.pt();   
    }
  );

  return jets;
}

std::vector<reco::CaloJet> Phase1L1TJetProducer::_buildJetsFromSeeds(const TH2F & caloGrid, const std::vector<std::tuple<int, int>> & seeds)
{
  // For each seed take a grid centered on the seed of the size specified by the user
  // Sum the pf in the grid, that will be the pt of the l1t jet. Eta and phi of the jet is taken from the seed.
  std::vector<reco::CaloJet> jets;
  for (const auto& seed: seeds)
  {
    reco::CaloJet jet = this -> _buildJetFromSeed(caloGrid, seed);
    if ( _debug ) std::cout << "jet pt-eta-phi: " << (float) jet.pt() << "\t" <<(float) jet.eta() << "\t" << (float) jet.phi() << std::endl;
    jets.push_back(jet);
  }

  // sort by pt
  std::sort(jets.begin(), jets.end(), 
    [](const reco::CaloJet& jet1, const reco::CaloJet& jet2) 
    {
      return jet1.pt() > jet2.pt();   
    }
  );
  return jets;
}


template <class Container>
void Phase1L1TJetProducer::_fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives, const unsigned int pfRegionIndex)
{
  //Filling the calo grid with the primitives
  // for ( const auto& primitiveIterator = triggerPrimitives.begin(); primitiveIterator != triggerPrimitives.end(); primitiveIterator++)
  for ( const auto primitiveIterator : triggerPrimitives )
    {
      // if(primitiveIterator->eta() >= 0 && primitiveIterator->eta() < 1.5 && primitiveIterator->phi() >= 0 && primitiveIterator->phi() < 0.7)
      // {
      //   std::cout << "input pt-eta-phi: " << (float) primitiveIterator -> pt() << "\t" <<(float) primitiveIterator -> eta() << "\t" << (float) primitiveIterator -> phi() << std::endl;
      // }

      if( primitiveIterator->phi() >= _phiLow &&
          primitiveIterator->phi() < _phiUp &&
          primitiveIterator->eta() >= _etaBinning.front() &&
          primitiveIterator->eta() < _etaBinning.back()){


        float eta = primitiveIterator->eta();
        float phi = primitiveIterator->phi();

        if ( _debug ) std::cout << "====== New input in _fillCaloGrid ======" << std::endl;
        if ( _debug )  std::cout << "Input pt, eta, phi : " << primitiveIterator->pt() << " " << eta << " " << phi << " " << caloGrid.FindBin( eta, phi) << std::endl;

        // Convert eta and phi to digitised quantities
        // Use to check behaviour with bin edges/pf region edges
        std::pair< double, double > pfRegionEdges = pfRegionEtaPhiLowEdges( pfRegionIndex );
        int digitisedEta = floor( ( eta - pfRegionEdges.second ) / _etalsb );
        int digitisedPhi = floor( phi / _philsb );

        if ( _debug )  std::cout << "Digi (int) eta, phi : " << digitisedEta << " " << digitisedPhi << std::endl;

        // If eta or phi is on a bin edge
        // Put in bin above, to match behaviour of HLS
        // Unless it's on the last bin of this pf region
        TAxis* etaAxis = caloGrid.GetXaxis();
        std::pair< double, double > pfRegionUpEdges = pfRegionEtaPhiUpEdges( pfRegionIndex );
        int digiEtaEdgeLastBinUp = floor( ( pfRegionUpEdges.second - pfRegionEdges.second ) / _etalsb );
        if ( digitisedEta >= digiEtaEdgeLastBinUp ) {
          digitisedEta = digiEtaEdgeLastBinUp-1;
        }
        else {
          for ( int i = 0; i < etaAxis->GetNbins(); ++i ) {
            int digiEdgeBinUp =  floor( ( etaAxis->GetBinUpEdge(i) - pfRegionEdges.second ) / _etalsb );
            if ( digiEdgeBinUp == digitisedEta ){
              digitisedEta += 1;
              if ( _debug ) std::cout << "Changed digi eta to : " << digitisedEta << std::endl;
            }
          }          
        }

        TAxis* phiAxis = caloGrid.GetYaxis();
        int digiPhiEdgeLastBinUp = floor( ( pfRegionUpEdges.first ) / _etalsb );
        if ( digitisedPhi >= digiPhiEdgeLastBinUp ) {
          digitisedPhi = digiPhiEdgeLastBinUp-1;
        }
        else {
          for ( int i = 0; i < phiAxis->GetNbins(); ++i ) {
            int digiEdgeBinUp =  floor( phiAxis->GetBinUpEdge(i) / _philsb );
            if ( digiEdgeBinUp == digitisedPhi ){
              digitisedPhi += 1;
              if ( _debug ) std::cout << "Changed digi phi to : " << digitisedPhi << std::endl;
            }
          }          
        }

        // Convert digitised eta and phi back to floatin point quantities with reduced precision
        eta = ( digitisedEta + 0.5 ) * _etalsb + pfRegionEdges.second;
        phi = ( digitisedPhi + 0.5 ) * _philsb;


        if ( _debug )  std::cout << "Digitised eta phi : " << eta << " " << phi << " " << caloGrid.FindBin( eta, phi) << std::endl;

        caloGrid.Fill( eta, phi, (float) primitiveIterator -> pt());
      }
    }
  return;
}

template <class Handle >
std::vector< std::vector< edm::Ptr< reco::Candidate > > > Phase1L1TJetProducer::_prepareInputsIntoRegions( const Handle triggerPrimitives )
{

  std::vector< std::vector< reco::CandidatePtr > > inputsInRegions{ _etaPFRegionEdges.size() * ( _phiPFRegionEdges.size() - 1 ) };

  for ( unsigned int i = 0; i < triggerPrimitives->size(); ++i ) {

    reco::CandidatePtr tp( triggerPrimitives, i );

    // Which phi region does this tp belong to
    auto it_phi = _phiPFRegionEdges.begin();
    if ( tp->phi() > *_phiPFRegionEdges.begin() ) {
      it_phi = std::upper_bound (_phiPFRegionEdges.begin(), _phiPFRegionEdges.end(), tp->phi()) - 1;
    }

    // Which eta region does this tp belong to
    auto it_eta = _etaPFRegionEdges.begin();
    if ( tp->eta() > *_etaPFRegionEdges.begin() ) {
      it_eta = std::upper_bound (_etaPFRegionEdges.begin(), _etaPFRegionEdges.end(), tp->eta()) - 1;
    }

    if ( it_phi != _phiPFRegionEdges.end() && it_eta != _etaPFRegionEdges.end() ) {
      auto phiRegion = it_phi - _phiPFRegionEdges.begin();
      auto etaRegion = it_eta - _etaPFRegionEdges.begin();

      inputsInRegions[ pfRegionIndex( phiRegion, etaRegion ) ].emplace_back( tp );
    }
  }

  // Truncate inputs in each pf region
  for ( auto& inputs : inputsInRegions ) {
    if ( inputs.size() > _maxInputsPerPFRegion ) {
      inputs.resize( _maxInputsPerPFRegion );
    }
  }

  return inputsInRegions;
}

unsigned int Phase1L1TJetProducer::pfRegionIndex( const unsigned int phiRegion, const unsigned int etaRegion ) {
  return etaRegion * ( _phiPFRegionEdges.size() - 1 ) + phiRegion;
}

std::pair< double, double> Phase1L1TJetProducer::pfRegionEtaPhiLowEdges( const unsigned int pfRegionIndex ) {
  unsigned int phiRegion = pfRegionIndex % ( _phiPFRegionEdges.size() - 1 );
  unsigned int etaRegion = ( pfRegionIndex - phiRegion ) / ( _phiPFRegionEdges.size() - 1 );
  return std::pair< double, double > { _phiPFRegionEdges.at( phiRegion ), _etaPFRegionEdges.at( etaRegion ) };
}

std::pair< double, double> Phase1L1TJetProducer::pfRegionEtaPhiUpEdges( const unsigned int pfRegionIndex ) {
  unsigned int phiRegion = pfRegionIndex % ( _phiPFRegionEdges.size() - 1 );
  unsigned int etaRegion = ( pfRegionIndex - phiRegion ) / ( _phiPFRegionEdges.size() - 1 );
  if ( phiRegion == _phiPFRegionEdges.size()-1 ) {
    return std::pair< double, double > { _phiPFRegionEdges.at( phiRegion ), _etaPFRegionEdges.at( etaRegion + 1 ) };
  }
  else if ( etaRegion == _etaPFRegionEdges.size()-1 ) {
    return std::pair< double, double > { _phiPFRegionEdges.at( phiRegion + 1), _etaPFRegionEdges.at( etaRegion ) };    
  }

  return std::pair< double, double > { _phiPFRegionEdges.at( phiRegion + 1 ), _etaPFRegionEdges.at( etaRegion + 1 ) };
}

// computes MET
// Takes grid used by jet finder and projects to 1D histogram of phi, bin contents are total pt in that phi bin
// the phi bin index is used to retrieve the sin-cos value from the LUT emulator
// the pt of the input is multiplied by that sin cos value to obtain px and py that is added to the total event px & py
// after all the inputs have been processed we compute the total pt of the event, and set that as MET
l1t::EtSum Phase1L1TJetProducer::_computeMET( const TH2F & caloGrid, double etaCut, l1t::EtSum::EtSumType sumType )
{
  const auto phiProjection = caloGrid.ProjectionY();

  unsigned int totalDigiPx{0};
  unsigned int totalDigiPy{0};

  for ( int i = 1; i < phiProjection->GetNbinsX()+1; ++i ) {
    double pt = phiProjection->GetBinContent( i );
    totalDigiPx += floor( floor( pt / 0.25 ) * _cosPhi[i-1] );
    totalDigiPy += floor( floor( pt / 0.25 ) * _sinPhi[i-1] );

    if ( _debug ) std::cout << i << "\t" << phiProjection->GetBinContent( i ) << "\t" << _sinPhi[i-1] << "\t" << floor( floor( pt / 0.25 ) * _sinPhi[i-1] ) << "\t" << _cosPhi[i-1] << "\t" << floor( floor( pt / 0.25 ) * _cosPhi[i-1] ) << std::endl;
  }

  double lMET = floor( sqrt(totalDigiPx * totalDigiPx + totalDigiPy * totalDigiPy) ) * 0.25;

  if ( _debug ) {
    std::cout << "Total digi px, py : " << totalDigiPx << " " << totalDigiPy << std::endl;
    std::cout << "Squared digi px, py : " << totalDigiPx * totalDigiPx << " " << totalDigiPy * totalDigiPy << std::endl;
    std::cout << "MET : " << floor( sqrt(totalDigiPx * totalDigiPx + totalDigiPy * totalDigiPy) ) << " " << lMET << std::endl;
  }

  // math::XYZTLorentzVector lMETVector( ( totalDigiPx ) * 0.25,  ( totalDigiPy ) * 0.25, 0, lMET);
  math::PtEtaPhiMLorentzVector lMETVector( lMET, 0, acos(totalDigiPx / ( lMET / 0.25 )), 0 );
  l1t::EtSum lMETSum(lMETVector, l1t::EtSum::EtSumType::kMissingEt, 0, 0, 0, 0 );
  // std::cout << lMETVector.pt() << " == " << lMET << "?" << std::endl;

  return lMETSum;
}

bool Phase1L1TJetProducer::_trimTower( const int etaIndex, const int phiIndex )
{

  int etaHalfSize = (int) this -> _jetIEtaSize/2;
  int phiHalfSize = (int) this -> _jetIPhiSize/2;

  if ( etaIndex == -etaHalfSize || etaIndex == etaHalfSize ) {
    if ( phiIndex <= -phiHalfSize+1 || phiIndex >= phiHalfSize-1 ) {
      return true;
    }
  }
  else if ( etaIndex == -etaHalfSize+1 || etaIndex == etaHalfSize-1 ) {
    if ( phiIndex == -phiHalfSize || phiIndex == phiHalfSize ) {
      return true;
    }    
  }

  return false;
}

void
Phase1L1TJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Phase1L1TJetProducer);
