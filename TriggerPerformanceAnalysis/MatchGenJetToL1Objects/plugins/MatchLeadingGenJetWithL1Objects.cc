// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchLeadingGenJetWithL1Objects
// Class:      MatchLeadingGenJetWithL1Objects
// 
/**\class MatchLeadingGenJetWithL1Objects MatchLeadingGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchLeadingGenJetWithL1Objects/plugins/MatchLeadingGenJetWithL1Objects.cc

 Description: Matches a gen jet with a L1T object and creates a tree

 Implementation:
    Receives tags of the various L1T objects
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 14:21:56 GMT
//
//


// system include files
#include <memory>
#include <csignal>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <utility>

#include "TriggerPerformanceAnalysis/MatchGenJetToL1Objects/interface/MatchingAlgorithms.hxx"

struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

struct TriggerObject {
  unsigned int id;
  float pt, eta, phi;
  int hwQual;
};

class MatchLeadingGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchLeadingGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchLeadingGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector< std::tuple < TriggerObject, Particle, float, int> >
    _matchGenJetWithL1Object(
      const reco::GenJet &,
      const edm::Handle<BXVector<T>> &,
      float = 0.25,
      bool = true
    );

    template <class TParticle, class TTrigger> // <3
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector<std::tuple<const TTrigger*, const TParticle*, float, int> > &
    );

    edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genParticleCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    edm::EDGetTokenT< std::vector< reco::CaloJet > > *_l1tJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Tau > > *_l1tTauCollectionTag;
    //!TODO: It should splitted in two different objects in the Stage-2
    edm::EDGetTokenT< BXVector< l1t::EGamma > > *_l1tEGammaCollectionTag;
    //!TODO: Like this
    //edm::EDGetTokenT< BXVector< l1t::Electron > > _l1tElectronCollectionTag;
    //edm::EDGetTokenT< BXVector< l1t::Gamma > > _l1tGammaCollectionTag;

    TTree * _l1tMuonLeadingGenJetTree;
    TTree * _l1tJetLeadingGenJetTree;
    TTree * _l1tEGammaLeadingGenJetTree;
    TTree * _l1tTauLeadingGenJetTree;
    TTree * _leadingGenJetTree;
    TTree * _leadingGenJetInFourPiTree;

    Particle _genJetParticle;
    TriggerObject _l1tObjectParticle;
    int _hwQual;
    float _deltaR2;
    int _matchingQuality;
};

MatchLeadingGenJetWithL1Objects::MatchLeadingGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TMuonTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tJetLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tEGammaLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TEGammaTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tTauLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TTauTree", "TTree with generator-level jet / L1T Tau information");
  this -> _leadingGenJetTree = fs -> make<TTree>("leadingGenJetTree", "TTree with generator-level jet information in detector space");
  this -> _leadingGenJetInFourPiTree = fs -> make<TTree>("_leadingGenJetInFourPiTree", "TTree with generator-level jet information in four pi");

  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  this -> _l1tJetLeadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  this -> _l1tEGammaLeadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_id", &(this -> _l1tObjectParticle.id), "l1tEGamma_id/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectParticle.pt), "l1tEGamma_pt/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectParticle.eta), "l1tEGamma_eta/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectParticle.phi), "l1tEGamma_phi/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_qual", &(this -> _l1tObjectParticle.hwQual), "l1tEGamma_qual/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  this -> _l1tTauLeadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_id", &(this -> _l1tObjectParticle.id), "l1tTau_id/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_pt", &(this -> _l1tObjectParticle.pt), "l1tTau_pt/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_eta", &(this -> _l1tObjectParticle.eta), "l1tTau_eta/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_phi", &(this -> _l1tObjectParticle.phi), "l1tTau_phi/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_qual", &(this -> _l1tObjectParticle.hwQual), "l1tTau_qual/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _leadingGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _leadingGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _leadingGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _leadingGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  
  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _leadingGenJetInFourPiTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _leadingGenJetInFourPiTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _leadingGenJetInFourPiTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _leadingGenJetInFourPiTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

}

void MatchLeadingGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genParticleCollectionTag not found." << std::endl;
    this -> _genParticleCollectionTag = NULL;
  }

  try
  {
    this -> _l1tMuonCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Muon > >(consumes< BXVector< l1t::Muon > > (iConfig.getParameter< edm::InputTag >("l1tMuonCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMuonCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tMuonCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::CaloJet > >(consumes< std::vector< reco::CaloJet > > (iConfig.getParameter< edm::InputTag >("l1tJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tJetCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tJetCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tEGammaCollectionTag = new edm::EDGetTokenT< BXVector < l1t::EGamma > >(consumes< BXVector< l1t::EGamma > > (iConfig.getParameter< edm::InputTag >("l1tEGammaCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tEGammaCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tEGammaCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tTauCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Tau > >(consumes< BXVector< l1t::Tau > > (iConfig.getParameter< edm::InputTag >("l1tTauCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tTauCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tTauCollectionTag = NULL;
  }
  
  return;
}

void MatchLeadingGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  if (this -> _l1tJetCollectionTag) delete this -> _l1tJetCollectionTag;
  if (this -> _l1tEGammaCollectionTag) delete this -> _l1tEGammaCollectionTag;
  if (this -> _l1tTauCollectionTag) delete this -> _l1tTauCollectionTag;
}

MatchLeadingGenJetWithL1Objects::~MatchLeadingGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchLeadingGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  //Looking for the leading gen jet

  float maxPt = 0;
  bool save = false;
  const reco::GenJet* leadingGenJetInBarrel = NULL;
  const reco::GenJet* leadingGenJet = NULL;

  for (auto genJetIterator = genJetCollectionHandle->begin(); genJetIterator != genJetCollectionHandle->end(); genJetIterator++)
  {
    // Only barrel leading l1t jets
    if ((genJetIterator->pt() > maxPt) && (genJetIterator->eta() < 5) && (genJetIterator->eta() > -5))
    {
      save = true;
      maxPt = genJetIterator->pt();
      leadingGenJetInBarrel = &(*genJetIterator);
      this->_genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
      this->_genJetParticle.pt = genJetIterator->pt();
      this->_genJetParticle.eta = genJetIterator->eta();
      this->_genJetParticle.phi = genJetIterator->phi();
    }
  }
  if (save)
    this->_leadingGenJetTree->Fill();

  maxPt = 0;
  save = false;

  for (auto genJetIterator = genJetCollectionHandle->begin(); genJetIterator != genJetCollectionHandle->end(); genJetIterator++)
  {
  // Only leading l1t jets
    if ((genJetIterator->pt() > maxPt))
    {
      save = true;
      maxPt = genJetIterator->pt();
      leadingGenJet = &(*genJetIterator);
      this->_genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
      this->_genJetParticle.pt = genJetIterator->pt();
      this->_genJetParticle.eta = genJetIterator->eta();
      this->_genJetParticle.phi = genJetIterator->phi();
    }
  }

  if (save)
    this->_leadingGenJetInFourPiTree->Fill();

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes

  if (leadingGenJet != NULL){
    if (this -> _l1tMuonCollectionTag)
    {
      // Muon veto. For muon matching
      edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
      iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);
      bool hasMuons = false;
      for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
      {
        if (abs(genParticleIterator->pdgId()) == 13)
          hasMuons = true;
      }

      if (!hasMuons)
      {
        //We look for the leading gen jet in four pi, the leading l1t muon in detector and we put them together
        edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
        float l1tMuonPtMax = 0;
        iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
        const l1t::Muon* leadingL1TMuon = NULL;
        for (auto l1tMuonIterator = l1tMuonCollectionHandle->begin(); l1tMuonIterator != l1tMuonCollectionHandle->end(); l1tMuonIterator++)
        {
          // leading muon in full detector space
          if (l1tMuonIterator->pt() > l1tMuonPtMax)
          {
            l1tMuonPtMax = l1tMuonIterator->pt();
            leadingL1TMuon = &(*l1tMuonIterator);
            
          }
        }

        // here we have leadingGenJet and leadingGenMuon
        if (leadingL1TMuon) {
          float dr2 = reco::deltaR2(*leadingGenJet, *leadingL1TMuon);

          std::tuple<const l1t::Muon *, const reco::GenJet *, float, int> l1tMuonGenJetPair;
          std::get<0>(l1tMuonGenJetPair) = leadingL1TMuon;
          std::get<1>(l1tMuonGenJetPair) = leadingGenJet;
          std::get<2>(l1tMuonGenJetPair) = dr2;
          std::get<3>(l1tMuonGenJetPair) = 0;

          std::vector<std::tuple<const l1t::Muon*,const  reco::GenJet*, float, int> > l1tMuonGenJetPairs;
          // if we have found a match let's add it.
          l1tMuonGenJetPairs.push_back(l1tMuonGenJetPair);
          this -> _fillTreeWithMatchedPairs(*(this->_l1tMuonLeadingGenJetTree), l1tMuonGenJetPairs);
        }
      }
    }
    
  }
  
  if (leadingGenJetInBarrel){

    if (this -> _l1tJetCollectionTag)
    {
      edm::Handle < std::vector< reco::CaloJet > > l1tJetCollectionHandle;
      iEvent.getByToken(*(this -> _l1tJetCollectionTag), l1tJetCollectionHandle);
      std::tuple<const reco::CaloJet *, const reco::GenJet *, float, int> l1tJetGenJetPair = MatchingAlgorithms::matchParticleWithObject<>(
          *leadingGenJetInBarrel,
          genJetCollectionHandle,
          l1tJetCollectionHandle,
          1,
          true);
      std::vector<std::tuple<const reco::CaloJet *, const reco::GenJet *, float, int>> l1tJetGenJetPairs;
      if (std::get<0>(l1tJetGenJetPair) != NULL) l1tJetGenJetPairs.push_back(l1tJetGenJetPair);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tJetLeadingGenJetTree), l1tJetGenJetPairs);
    }

    if (this -> _l1tEGammaCollectionTag)
    {
      edm::Handle < BXVector< l1t::EGamma > > l1tEGammaCollectionHandle;
      iEvent.getByToken(*(this -> _l1tEGammaCollectionTag), l1tEGammaCollectionHandle);
      std::tuple<const l1t::EGamma *, const reco::GenJet *, float, int> l1tEGammaGenJetPair = MatchingAlgorithms::matchParticleWithL1Object<>(
          *leadingGenJetInBarrel,
          genJetCollectionHandle,
          l1tEGammaCollectionHandle,
          1,
          true);
      std::vector<std::tuple<const l1t::EGamma *, const reco::GenJet *, float, int>> l1tEGammaGenJetPairs;
      if (std::get<0>(l1tEGammaGenJetPair) != NULL) l1tEGammaGenJetPairs.push_back(l1tEGammaGenJetPair);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tEGammaLeadingGenJetTree), l1tEGammaGenJetPairs);
    }
    
    if (this -> _l1tTauCollectionTag)
    {
      edm::Handle < BXVector< l1t::Tau > > l1tTauCollectionHandle;
      iEvent.getByToken(*(this -> _l1tTauCollectionTag), l1tTauCollectionHandle);
      std::tuple<const l1t::Tau *, const reco::GenJet *, float, int> l1tTauGenJetPair = MatchingAlgorithms::matchParticleWithL1Object<>(
          *leadingGenJetInBarrel,
          genJetCollectionHandle,
          l1tTauCollectionHandle,
          1,
          true);
      std::vector<std::tuple<const l1t::Tau *, const reco::GenJet *, float, int>> l1tTauGenJetPairs;
      if (std::get<0>(l1tTauGenJetPair) != NULL) l1tTauGenJetPairs.push_back(l1tTauGenJetPair);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tTauLeadingGenJetTree), l1tTauGenJetPairs);
    }
  }
}

template <class TParticle, class TTrigger> // <3
void MatchLeadingGenJetWithL1Objects::_fillTreeWithMatchedPairs(
    TTree &aTree,
    const std::vector< std::tuple<const TTrigger*,const  TParticle*, float, int>> &matchedL1TObjectJetPairs)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const TTrigger & l1tObject = *(std::get<0>(matchTuple));
    const TParticle & genJet = *(std::get<1>(matchTuple));
    float deltaR2 = std::get<2>(matchTuple);
    int matchingQuality = std::get<3>(matchTuple);
    this -> _genJetParticle.id = 0;
    this -> _genJetParticle.pt = genJet.pt();
    this -> _genJetParticle.eta = genJet.eta();
    this -> _genJetParticle.phi = genJet.phi();
    this -> _l1tObjectParticle.id = 0;
    this -> _l1tObjectParticle.pt = l1tObject.pt();
    this -> _l1tObjectParticle.eta = l1tObject.eta();
    this -> _l1tObjectParticle.phi = l1tObject.phi();
    this -> _l1tObjectParticle.hwQual = 0;
    this -> _deltaR2 = deltaR2;
    this -> _matchingQuality = matchingQuality;
    aTree.Fill();
  }
}

// Considers a genJet and looks for the closest l1tObject within that radius (1 genjet can have only one l1tobject)

template <class T> // <3
const std::vector <std::tuple < TriggerObject, Particle, float, int > > 
MatchLeadingGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const reco::GenJet & genJet,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle,
  float dr2Min,
  bool crossMatch
)
{

  std::vector< std::tuple < TriggerObject, Particle , float, int > > l1tObjectGenJetPairs;
  // for each object in the genJet collection we look for the closest l1tobject in a wide range
  
  bool foundMatch = false;
  std::tuple<TriggerObject, Particle, float, int > l1tObjectGenJetPair;
  
  Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
  matchedGenJet.id = 0;
  matchedGenJet.pt = genJet.pt();
  matchedGenJet.phi = genJet.phi();
  matchedGenJet.eta = genJet.eta();

  for (
        typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin(0);
        bx0Iterator != l1tObjectCollectionHandle->end(0);
        bx0Iterator++
  )
  {
    float dr2 = reco::deltaR2(*bx0Iterator, genJet);

    if ((dr2 < dr2Min))
    {
      TriggerObject &matchedL1TObject = std::get<0>(l1tObjectGenJetPair);
      matchedL1TObject.id = (bx0Iterator - l1tObjectCollectionHandle->begin(0));
      matchedL1TObject.pt = bx0Iterator->pt();
      matchedL1TObject.phi = bx0Iterator->phi();
      matchedL1TObject.eta = bx0Iterator->eta();
      matchedL1TObject.hwQual = bx0Iterator->hwQual();

      std::get<2>(l1tObjectGenJetPair) = dr2;
      dr2Min = dr2;
      foundMatch = true;
    }
  
  } 
  //if we have found a compatible l1t object we push the object-jet pair in a vector which will be our result
  if (foundMatch)
    l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchLeadingGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchLeadingGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchLeadingGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchLeadingGenJetWithL1Objects);
