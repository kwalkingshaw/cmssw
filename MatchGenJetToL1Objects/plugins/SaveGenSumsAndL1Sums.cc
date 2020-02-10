// -*- C++ -*-
//
// Package:    caloJetConvolutionCurves/SaveGenSumsAndL1Sums
// Class:      SaveGenSumsAndL1Sums
// 
/**\class SaveGenSumsAndL1Sums SaveGenSumsAndL1Sums.cc caloJetConvolutionCurves/SaveGenSumsAndL1Sums/plugins/SaveGenSumsAndL1Sums.cc
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
#include <csignal>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <utility>


// struct Particle {
//   unsigned int id;
//   float pt, eta, phi;
// };

// struct CaloJet {
//   unsigned int id;
//   float pt, eta, phi;
//   float pileup;
// };

class SaveGenSumsAndL1Sums : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit SaveGenSumsAndL1Sums(const edm::ParameterSet&);
    ~SaveGenSumsAndL1Sums();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();
    double _computeHT(const std::vector<reco::GenJet>& jetVector);

    edm::EDGetTokenT< std::vector< reco::GenMET > > *_genMETCollectionTag;
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tMETCollectionTag;
    
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tHTCollectionTag;
    
    TTree * _genMETL1TMETTree;
    TTree * _genHTL1THTTree;

    float _genMET, _l1tMET, _genHT, _l1tHT;
    
    double _htPtThreshold;

};

SaveGenSumsAndL1Sums::SaveGenSumsAndL1Sums(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  if ((this -> _genMETCollectionTag) && (this -> _l1tMETCollectionTag))
  {
    this -> _genMETL1TMETTree = fs -> make<TTree>("genMETL1TMETTree", "TTree with generator-level / L1T MET information");
    this -> _genMETL1TMETTree -> Branch("genMET", &(this -> _genMET), "genMET/F");
    this -> _genMETL1TMETTree -> Branch("l1tMET", &(this -> _l1tMET), "l1TMET/F");
  }
  
  if ((this -> _genJetCollectionTag) && (this -> _l1tHTCollectionTag))
  {
    this -> _genHTL1THTTree = fs -> make<TTree>("genHTL1THTTree", "TTree with generator-level / L1T HT information");
    this -> _genHTL1THTTree -> Branch("genHT", &(this -> _genHT), "genHT/F");
    this -> _genHTL1THTTree -> Branch("l1tHT", &(this -> _l1tHT), "l1tHT/F");
  }
}

void SaveGenSumsAndL1Sums::_getTokens(const edm::ParameterSet& iConfig)
{
  
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genMETCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenMET > >(consumes< std::vector< reco::GenMET > > (iConfig.getParameter< edm::InputTag >("genMETCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> MET configuration not found, proceeding without adding MET info to tree" << std::endl;
    this -> _genMETCollectionTag = NULL;
  }

  try
  {
    this -> _l1tMETCollectionTag = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tMETCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> MET configuration not found, proceeding without adding MET info to tree" << std::endl;
    this -> _l1tMETCollectionTag = NULL;
  }

  try
  {
    this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> HT configuration not found, proceeding without adding HT info to tree" << std::endl;
    this -> _genJetCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tHTCollectionTag = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tHTCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> HT configuration not found, proceeding without adding HT info to tree" << std::endl;
    this -> _l1tHTCollectionTag = NULL;
  }

  try
  {
    this -> _htPtThreshold = iConfig.getParameter<double>("htPtThreshold");
  } catch (std::exception const & ex) 
  {
    this -> _htPtThreshold = 30;
  }

  return;
}

void SaveGenSumsAndL1Sums::_freeTokens()
{
  if (this -> _genMETCollectionTag) delete this -> _genMETCollectionTag;
  if (this -> _l1tMETCollectionTag) delete this -> _l1tMETCollectionTag;
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _l1tHTCollectionTag) delete this -> _l1tHTCollectionTag;
  return;
}

SaveGenSumsAndL1Sums::~SaveGenSumsAndL1Sums()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SaveGenSumsAndL1Sums::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  
  if ((this -> _genMETCollectionTag) && (this -> _l1tMETCollectionTag))
  {
    edm::Handle < std::vector< reco::GenMET> > lGenMETCollectionHandle;
    edm::Handle < BXVector<l1t::EtSum> > lL1TMETCollectionHandle;
    iEvent.getByToken(*(this -> _genMETCollectionTag), lGenMETCollectionHandle);
    iEvent.getByToken(*(this -> _l1tMETCollectionTag), lL1TMETCollectionHandle);
    //access gen and l1t met data
    this -> _genMET = lGenMETCollectionHandle -> front().pt();
    const BXVector<l1t::EtSum> & lL1TMETCollection = *lL1TMETCollectionHandle;
    for (const l1t::EtSum & lL1TMET: lL1TMETCollection) if (lL1TMET.getType() == l1t::EtSum::EtSumType::kMissingEt) this -> _l1tMET = lL1TMET.pt();
    this -> _genMETL1TMETTree -> Fill();
  }

  if ((this -> _genJetCollectionTag) && (this -> _l1tHTCollectionTag))
  {
    edm::Handle < std::vector< reco::GenJet> > lGenJetCollectionHandle;
    edm::Handle < BXVector<l1t::EtSum> > lL1THTCollectionHandle;
    iEvent.getByToken(*(this -> _genJetCollectionTag), lGenJetCollectionHandle);
    iEvent.getByToken(*(this -> _l1tHTCollectionTag), lL1THTCollectionHandle);
    // access l1t HT data 
    const BXVector<l1t::EtSum> & lL1THTCollection = *lL1THTCollectionHandle;
    for (const l1t::EtSum & lL1THT: lL1THTCollection) if (lL1THT.getType() == l1t::EtSum::EtSumType::kTotalHt) this -> _l1tHT = lL1THT.pt();
    // compute gen HT
    this -> _genHT = this -> _computeHT(*lGenJetCollectionHandle);
    this -> _genHTL1THTTree -> Fill();
  }

}

double SaveGenSumsAndL1Sums::_computeHT(const std::vector<reco::GenJet>& jetVector) 
{
  double lHT = 0;

  for (const auto & jet: jetVector)
  {
    lHT += (jet.pt() >= this -> _htPtThreshold) ? jet.pt() : 0;
  }

  return lHT;
}

// ------------ method called once each job just before starting event loop  ------------
void SaveGenSumsAndL1Sums::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SaveGenSumsAndL1Sums::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SaveGenSumsAndL1Sums::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SaveGenSumsAndL1Sums);