// -*- C++ -*-
//
// Package:    caloJetConvolutionCurves/MatchGenJetToRecoCaloJet
// Class:      MatchGenJetToRecoCaloJet
// 
/**\class MatchGenJetToRecoCaloJet MatchGenJetToRecoCaloJet.cc caloJetConvolutionCurves/MatchGenJetToRecoCaloJet/plugins/MatchGenJetToRecoCaloJet.cc
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
#include <signal.h>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/JetReco/interface/CaloJet.h"

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

struct CaloJet {
  unsigned int id;
  float pt, eta, phi;
  float pileup;
};

class MatchGenJetToRecoCaloJet : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchGenJetToRecoCaloJet(const edm::ParameterSet&);
    ~MatchGenJetToRecoCaloJet();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector< std::tuple < CaloJet, Particle, float> >
    _matchGenJetWithObject
    (
      const edm::Handle<std::vector<reco::GenJet>> &,
      const edm::Handle<std::vector<T>> &,
      float = 0.25
    );

    template <class TParticle, class TTrigger> // <3
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector < std::tuple < const TTrigger*, const TParticle*, float, int > > &
    );

    template <class TParticle, class TTrigger> // <3
    const std::vector< std::tuple<const TTrigger*, const TParticle*, float, int> >
    _matchParticleWithObject
    (
      const edm::Handle<std::vector<TParticle>> &,
      const edm::Handle<std::vector<TTrigger>> &,
      float = 0.25,
      bool = true
    );

    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< std::vector< reco::CaloJet > > *_caloJetCollectionTag;

    TTree * _caloJetGenJetTree;
    TTree * _caloJetTree;
    TTree * _caloLeadingJetTree;
    TTree * _genJetTree;

    Particle _genJetParticle;
    CaloJet _caloJetParticle;
    float _deltaR2;
    int _matchingQuality;
};

MatchGenJetToRecoCaloJet::MatchGenJetToRecoCaloJet(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _caloJetGenJetTree = fs -> make<TTree>("matchedCaloJetGenJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _caloJetTree = fs -> make<TTree>("caloJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _caloLeadingJetTree = fs -> make<TTree>("caloLeadingJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _genJetTree = fs -> make<TTree>("genJetTree", "TTree with generator-level jet information");

  this -> _caloJetGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _caloJetGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _caloJetGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _caloJetGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _caloJetGenJetTree -> Branch("caloJet_id", &(this -> _caloJetParticle.id), "caloJet_id/i");
  this -> _caloJetGenJetTree -> Branch("caloJet_pt", &(this -> _caloJetParticle.pt), "caloJet_pt/F");
  this -> _caloJetGenJetTree -> Branch("caloJet_eta", &(this -> _caloJetParticle.eta), "caloJet_eta/F");
  this -> _caloJetGenJetTree -> Branch("caloJet_phi", &(this -> _caloJetParticle.phi), "caloJet_phi/F");
  this -> _caloJetGenJetTree -> Branch("caloJet_pileup", &(this -> _caloJetParticle.pileup), "caloJet_pileup/F");
  this -> _caloJetGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");
  this -> _caloJetGenJetTree -> Branch("matchingQuality", &(this -> _matchingQuality), "matchingQuality/i");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _genJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _genJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _genJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _caloJetTree -> Branch("caloJet_id", &(this -> _caloJetParticle.id), "caloJet_id/i");
  this -> _caloJetTree -> Branch("caloJet_pt", &(this -> _caloJetParticle.pt), "caloJet_pt/F");
  this -> _caloJetTree -> Branch("caloJet_eta", &(this -> _caloJetParticle.eta), "caloJet_eta/F");
  this -> _caloJetTree -> Branch("caloJet_phi", &(this -> _caloJetParticle.phi), "caloJet_phi/F");
  this -> _caloJetTree -> Branch("caloJet_pileup", &(this -> _caloJetParticle.pileup), "caloJet_pileup/F");

  this -> _caloLeadingJetTree -> Branch("caloJet_id", &(this -> _caloJetParticle.id), "caloJet_id/i");
  this -> _caloLeadingJetTree -> Branch("caloJet_pt", &(this -> _caloJetParticle.pt), "caloJet_pt/F");
  this -> _caloLeadingJetTree -> Branch("caloJet_eta", &(this -> _caloJetParticle.eta), "caloJet_eta/F");
  this -> _caloLeadingJetTree -> Branch("caloJet_phi", &(this -> _caloJetParticle.phi), "caloJet_phi/F");
  this -> _caloLeadingJetTree -> Branch("caloJet_pileup", &(this -> _caloJetParticle.pileup), "caloJet_pileup/F");

}

void MatchGenJetToRecoCaloJet::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _caloJetCollectionTag = new edm::EDGetTokenT< std::vector < reco::CaloJet > >(consumes< std::vector< reco::CaloJet > > (iConfig.getParameter< edm::InputTag >("caloJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> caloJetCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _caloJetCollectionTag = NULL;
  }
  return;
}

void MatchGenJetToRecoCaloJet::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _caloJetCollectionTag) delete this -> _caloJetCollectionTag;
}

MatchGenJetToRecoCaloJet::~MatchGenJetToRecoCaloJet()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchGenJetToRecoCaloJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  
  if (this -> _caloJetCollectionTag)
  {
    edm::Handle < std::vector< reco::CaloJet > > caloJetCollectionHandle;
    iEvent.getByToken(*(this -> _caloJetCollectionTag), caloJetCollectionHandle);
    const auto caloJetGenJetPairs = this -> _matchParticleWithObject<>(genJetCollectionHandle, caloJetCollectionHandle, 1);
    this -> _fillTreeWithMatchedPairs(*(this -> _caloJetGenJetTree), caloJetGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto caloJetIterator = caloJetCollectionHandle -> begin(); caloJetIterator != caloJetCollectionHandle -> end(); caloJetIterator++ )
    {
      // Only barrel leading l1t jets
      if ((caloJetIterator->pt() > maxPt) && (caloJetIterator->eta() < 1.44) && (caloJetIterator->eta() > -1.44))
      {
        save = true;        
        maxPt = caloJetIterator -> pt();
        this -> _caloJetParticle.pileup = caloJetIterator -> pileup();      
        this -> _caloJetParticle.id = (caloJetIterator - caloJetCollectionHandle->begin());
        this -> _caloJetParticle.pt = caloJetIterator -> pt();
        this -> _caloJetParticle.eta = caloJetIterator -> eta();
        this -> _caloJetParticle.phi = caloJetIterator -> phi();
      }
    }
    if (save) this -> _caloLeadingJetTree -> Fill();
    
    for (auto caloJetIterator = caloJetCollectionHandle -> begin(); caloJetIterator != caloJetCollectionHandle -> end(); caloJetIterator++ )
    {
      this -> _caloJetParticle.id = (caloJetIterator - caloJetCollectionHandle->begin());
      this -> _caloJetParticle.pt = caloJetIterator -> pt();
      this -> _caloJetParticle.eta = caloJetIterator -> eta();
      this -> _caloJetParticle.phi = caloJetIterator -> phi();
      this -> _caloJetParticle.pileup = caloJetIterator -> pileup();            
      this -> _caloJetTree -> Fill();
      //std::cout << "jet pt-eta-phi " << this->_caloJetParticle.pt << "\t" << this->_caloJetParticle.eta << "\t" << this->_caloJetParticle.phi << std::endl;
    }
    
  }

  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    this -> _genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
    this -> _genJetParticle.pt = genJetIterator -> pt();
    //std::cout << "Jet pt " << genJetIterator -> pt() << std::endl;
    this -> _genJetParticle.eta = genJetIterator -> eta();
    this -> _genJetParticle.phi = genJetIterator -> phi();
    this -> _genJetTree -> Fill();
  }


}

template<class TParticle, class TTrigger>
void
MatchGenJetToRecoCaloJet::_fillTreeWithMatchedPairs
(
  TTree & aTree,
  const std::vector < std::tuple <const TTrigger*, const TParticle*, float, int > > & matchedL1TObjectJetPairs
)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const TTrigger & l1tObject = *(std::get<0>(matchTuple));
    const TParticle & genJet = *(std::get<1>(matchTuple));
    float deltaR2 = std::get<2>(matchTuple);
    float matchingQuality = std::get<3>(matchTuple);
    this -> _genJetParticle.id = 0;
    this -> _genJetParticle.pt = genJet.pt();
    this -> _genJetParticle.eta = genJet.eta();
    this -> _genJetParticle.phi = genJet.phi();
    this -> _caloJetParticle.id = 0;
    this -> _caloJetParticle.pt = l1tObject.pt();
    this -> _caloJetParticle.eta = l1tObject.eta();
    this -> _caloJetParticle.phi = l1tObject.phi();
    this -> _caloJetParticle.pileup = l1tObject.pileup();
    this -> _deltaR2 = deltaR2;
    this -> _matchingQuality = matchingQuality;
    
    aTree.Fill();
    
  }
}

template <class TParticle, class TTrigger> // <3
const std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > 
MatchGenJetToRecoCaloJet::_matchParticleWithObject
(
  const edm::Handle<std::vector<TParticle>>& particleCollectionHandle,
  const edm::Handle<std::vector<TTrigger>>& l1tObjectCollectionHandle,
  float dr2Min,
  bool crossMatch
)
{
  std::vector <std::tuple<const TTrigger*, const TParticle*, float, int> > l1tObjectParticlePairs;
  // for each object in the particle collection we look for the closest l1tobject in a wide range
  for (auto particleIterator = particleCollectionHandle -> begin(); particleIterator != particleCollectionHandle -> end(); particleIterator++ )
  {
    std::tuple<const TTrigger *, const TParticle *, float, int> l1tObjectParticlePair = 
      MatchingAlgorithms::matchParticleWithObject<>
      (
        *particleIterator,
        particleCollectionHandle,
        l1tObjectCollectionHandle,
        dr2Min,
        crossMatch
      );
    // if we have found a match let's add it.
    if (std::get<0>(l1tObjectParticlePair) != NULL) l1tObjectParticlePairs.push_back(l1tObjectParticlePair);
  }
  // return the matched stuff
  return l1tObjectParticlePairs;
}

// Considers a genJet and looks for the closest l1tObject within that radius (1 genjet can have only one l1tobject)

template <class T> // <3
const std::vector <std::tuple < CaloJet, Particle, float > > 
MatchGenJetToRecoCaloJet::_matchGenJetWithObject
(
  const edm::Handle < std::vector< reco::GenJet > > & genJetCollectionHandle,
  const edm::Handle < std::vector < T > > & l1tObjectCollectionHandle,
  float dr2Min
)
{

  std::vector< std::tuple < CaloJet, Particle , float > > l1tObjectGenJetPairs;
  // for each object in the genJet collection we look for the closest l1tobject in a wide range
  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {    
    bool foundMatch = false;
    std::tuple<CaloJet, Particle, float > l1tObjectGenJetPair;
    
    Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
    matchedGenJet.id = (genJetIterator - genJetCollectionHandle -> begin());
    matchedGenJet.pt = genJetIterator -> pt();
    matchedGenJet.phi = genJetIterator -> phi();
    matchedGenJet.eta = genJetIterator -> eta();

    for (
          typename std::vector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin();
          bx0Iterator != l1tObjectCollectionHandle->end();
          bx0Iterator++
    )
    {
      float dr2 = reco::deltaR2(*bx0Iterator, *genJetIterator);

      if ((dr2 < dr2Min))
      {
        CaloJet &matchedL1TObject = std::get<0>(l1tObjectGenJetPair);
        matchedL1TObject.id = (bx0Iterator - l1tObjectCollectionHandle->begin());
        matchedL1TObject.pt = bx0Iterator->pt();
        matchedL1TObject.phi = bx0Iterator->phi();
        matchedL1TObject.eta = bx0Iterator->eta();

        std::get<2>(l1tObjectGenJetPair) = dr2;
        dr2Min = dr2;
        foundMatch = true;
      }
    
    } 
    //if we have found a compatible l1t object we push the object-jet pair in a vector which will be our result
    if (foundMatch)
      l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  }

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchGenJetToRecoCaloJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchGenJetToRecoCaloJet::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchGenJetToRecoCaloJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchGenJetToRecoCaloJet);