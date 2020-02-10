// -*- C++ -*-
//
// Package:    caloJetConvolutionCurves/SaveL1Sums
// Class:      SaveL1Sums
// 
/**\class SaveL1Sums SaveL1Sums.cc caloJetConvolutionCurves/SaveL1Sums/plugins/SaveL1Sums.cc
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

class SaveL1Sums : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit SaveL1Sums(const edm::ParameterSet&);
    ~SaveL1Sums();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tMETCollectionTag1;
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tMETCollectionTag2;
    
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tHTCollectionTag1;
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tHTCollectionTag2;
    
    TTree * _l1tMETsTree;
    TTree * _l1tHTsTree;

    float _l1tMET1, _l1tMET2, _l1tHT1, _l1tHT2;

    std::string _l1tMETBranchName1, _l1tMETBranchName2;
    std::string _l1tHTBranchName1, _l1tHTBranchName2;
    
};

SaveL1Sums::SaveL1Sums(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  if ((this -> _l1tMETCollectionTag1) && (this -> _l1tMETCollectionTag2))
  {
    std::ostringstream lOSSL1TMETBranchName1;
    lOSSL1TMETBranchName1 << this -> _l1tMETBranchName1 << "/F";
    std::ostringstream lOSSL1TMETBranchName2;
    lOSSL1TMETBranchName2 << this -> _l1tMETBranchName2 << "/F";

    this -> _l1tMETsTree = fs -> make<TTree>("l1tMETsTree", "TTree with L1T-L1T MET information");
    this -> _l1tMETsTree -> Branch(this -> _l1tMETBranchName1.c_str(), &(this -> _l1tMET1), lOSSL1TMETBranchName1.str().c_str());
    this -> _l1tMETsTree -> Branch(this -> _l1tMETBranchName2.c_str(), &(this -> _l1tMET2), lOSSL1TMETBranchName2.str().c_str());
  }
  
  if ((this -> _l1tHTCollectionTag1) && (this -> _l1tHTCollectionTag2))
  {
    std::ostringstream lOSSL1THTBranchName1;
    lOSSL1THTBranchName1 << this -> _l1tHTBranchName1 << "/F";
    std::ostringstream lOSSL1THTBranchName2;
    lOSSL1THTBranchName2 << this -> _l1tHTBranchName2 << "/F";

    this -> _l1tHTsTree = fs -> make<TTree>("l1tHTsTree", "TTree with L1T-L1T HT information");
    this -> _l1tHTsTree -> Branch(this -> _l1tHTBranchName1.c_str(), &(this -> _l1tHT1), lOSSL1THTBranchName1.str().c_str());
    this -> _l1tHTsTree -> Branch(this -> _l1tHTBranchName2.c_str(), &(this -> _l1tHT2), lOSSL1THTBranchName2.str().c_str());
  }
}

void SaveL1Sums::_getTokens(const edm::ParameterSet& iConfig)
{
  
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _l1tMETCollectionTag1 = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tMETCollectionTag1")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMETCollectionTag1 not found, proceeding without adding MET info to tree" << std::endl;
    this -> _l1tMETCollectionTag1 = NULL;
  }

  try
  {
    this -> _l1tMETCollectionTag2 = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tMETCollectionTag2")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMETCollectionTag2 not found, proceeding without adding MET info to tree" << std::endl;
    this -> _l1tMETCollectionTag2 = NULL;
  }

  try
  {
    this -> _l1tMETBranchName1 = iConfig.getParameter<std::string>("l1tMETBranchName1");
    this -> _l1tMETBranchName2 = iConfig.getParameter<std::string>("l1tMETBranchName2");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMETBranchName1 or l1tMETBranchName2 not found, proceeding without adding MET info to tree" << std::endl;
    if (this -> _l1tMETCollectionTag1) 
    {
      delete this -> _l1tMETCollectionTag1;
      this -> _l1tMETCollectionTag1 = NULL;
    }
    if (this -> _l1tMETCollectionTag2) 
    {
      delete this -> _l1tMETCollectionTag2;
      this -> _l1tMETCollectionTag2 = NULL;
    }
  }

  try
  {
    this -> _l1tHTCollectionTag1 = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tHTCollectionTag1")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tHTCollectionTag1 not found, proceeding without adding HT info to tree" << std::endl;
    this -> _l1tHTCollectionTag1 = NULL;
  }
  
  try
  {
    this -> _l1tHTCollectionTag2 = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tHTCollectionTag2")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tHTCollectionTag2 not found, proceeding without adding HT info to tree" << std::endl;
    this -> _l1tHTCollectionTag2 = NULL;
  }

  try
  {
    this -> _l1tHTBranchName1 = iConfig.getParameter<std::string>("l1tHTBranchName1");
    this -> _l1tHTBranchName2 = iConfig.getParameter<std::string>("l1tHTBranchName2");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tHTBranchName1 or l1tHTBranchName2 not found, proceeding without adding HT info to tree" << std::endl;
    if (this -> _l1tHTCollectionTag1) 
    {
      delete this -> _l1tHTCollectionTag1;
      this -> _l1tHTCollectionTag1 = NULL;
    }
    if (this -> _l1tHTCollectionTag2) 
    {
      delete this -> _l1tHTCollectionTag2;
      this -> _l1tHTCollectionTag2 = NULL;
    }
  }

  return;
}

void SaveL1Sums::_freeTokens()
{
  if (this -> _l1tMETCollectionTag1) delete this -> _l1tMETCollectionTag1;
  if (this -> _l1tMETCollectionTag2) delete this -> _l1tMETCollectionTag2;
  if (this -> _l1tHTCollectionTag1) delete this -> _l1tHTCollectionTag1;
  if (this -> _l1tHTCollectionTag2) delete this -> _l1tHTCollectionTag2;
  return;
}

SaveL1Sums::~SaveL1Sums()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SaveL1Sums::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  
  if ((this -> _l1tMETCollectionTag1) && (this -> _l1tMETCollectionTag2))
  {
    edm::Handle < BXVector<l1t::EtSum> > lL1TMETCollectionHandle1;
    edm::Handle < BXVector<l1t::EtSum> > lL1TMETCollectionHandle2;
    iEvent.getByToken(*(this -> _l1tMETCollectionTag1), lL1TMETCollectionHandle1);
    iEvent.getByToken(*(this -> _l1tMETCollectionTag2), lL1TMETCollectionHandle2);
    //access gen and l1t met data
    const BXVector<l1t::EtSum> & lL1TMETCollection1 = *lL1TMETCollectionHandle1;
    const BXVector<l1t::EtSum> & lL1TMETCollection2 = *lL1TMETCollectionHandle2;
    for (const l1t::EtSum & lL1TMET: lL1TMETCollection1) if (lL1TMET.getType() == l1t::EtSum::EtSumType::kMissingEt) this -> _l1tMET1 = lL1TMET.pt();
    for (const l1t::EtSum & lL1TMET: lL1TMETCollection2) if (lL1TMET.getType() == l1t::EtSum::EtSumType::kMissingEt) this -> _l1tMET2 = lL1TMET.pt();
    this -> _l1tMETsTree -> Fill();
  }

  if ((this -> _l1tHTCollectionTag1) && (this -> _l1tHTCollectionTag2))
  {
    edm::Handle < BXVector<l1t::EtSum> > lL1THTCollectionHandle1;
    edm::Handle < BXVector<l1t::EtSum> > lL1THTCollectionHandle2;
    iEvent.getByToken(*(this -> _l1tHTCollectionTag1), lL1THTCollectionHandle1);
    iEvent.getByToken(*(this -> _l1tHTCollectionTag2), lL1THTCollectionHandle2);
    // access l1t HT data 
    const BXVector<l1t::EtSum> & lL1THTCollection1 = *lL1THTCollectionHandle1;
    const BXVector<l1t::EtSum> & lL1THTCollection2 = *lL1THTCollectionHandle2;
    for (const l1t::EtSum & lL1THT: lL1THTCollection1) if (lL1THT.getType() == l1t::EtSum::EtSumType::kTotalHt) this -> _l1tHT1 = lL1THT.pt();
    for (const l1t::EtSum & lL1THT: lL1THTCollection2) if (lL1THT.getType() == l1t::EtSum::EtSumType::kTotalHt) this -> _l1tHT2 = lL1THT.pt();
    this -> _l1tHTsTree -> Fill();
  }

}

// ------------ method called once each job just before starting event loop  ------------
void SaveL1Sums::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SaveL1Sums::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SaveL1Sums::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SaveL1Sums);