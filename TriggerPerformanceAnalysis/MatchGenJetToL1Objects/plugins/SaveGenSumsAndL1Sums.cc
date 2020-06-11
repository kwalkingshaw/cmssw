// -*- C++ -*-
//
// Package:    caloJetConvolutionCurves/SaveGenSumsAndL1Sums
// Class:      SaveGenSumsAndL1Sums
// 
/**\class SaveGenSumsAndL1Sums SaveGenSumsAndL1Sums.cc caloJetConvolutionCurves/SaveGenSumsAndL1Sums/plugins/SaveGenSumsAndL1Sums.cc
 Description: Creates a tree that saves generator-level and trigger-level sums
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
    double _computeMHT(const std::vector<reco::GenJet> & jetVector);

    //tag of the generator-level met
    edm::EDGetTokenT< std::vector< reco::GenMET > > *_genMETCollectionTag;
    //tag of the trigger-level met
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tMETCollectionTag;
    //tag of the generator-level jets HT and MHT can be computed from
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    //tag of the trigger-level ht
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tHTCollectionTag;
    //tag of the trigger-level mht
    edm::EDGetTokenT< BXVector<l1t::EtSum> > *_l1tMHTCollectionTag;
    
    //TTree holding gen-l1t MET pairs
    TTree * _genMETL1TMETTree;
    //TTree holding gen-l1t HT pairs
    TTree * _genHTL1THTTree;
    //TTree holding gen-l1t MHT pairs
    TTree * _genMHTL1TMHTTree;
    

    //memory area where to store gen and l1t sums before saving them into the tree
    float _genMET, _l1tMET, _genHT, _l1tHT, _genMHT, _l1tMHT;
    
    double _htPtThreshold;

};

//Obtains the token to access the inputs, initialises the TTree based on inputs
SaveGenSumsAndL1Sums::SaveGenSumsAndL1Sums(const edm::ParameterSet& iConfig)
{  
  //obtaining the tokens to access the input data
  this -> _getTokens(iConfig);
  //accessing the file where the TTrees will be saved to
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  //Creating a TTree where to save gen-l1t MET, if I have received the MET tags
  if ((this -> _genMETCollectionTag) && (this -> _l1tMETCollectionTag))
  {
    this -> _genMETL1TMETTree = fs -> make<TTree>("genMETL1TMETTree", "TTree with generator-level / L1T MET information");
    this -> _genMETL1TMETTree -> Branch("genMET", &(this -> _genMET), "genMET/F");
    this -> _genMETL1TMETTree -> Branch("l1tMET", &(this -> _l1tMET), "l1TMET/F");
  }
  
  //Creating a TTree where to save gen-l1t HT, if I have received the HT tag and the gen-jet tag required to compute it
  if ((this -> _genJetCollectionTag) && (this -> _l1tHTCollectionTag))
  {
    this -> _genHTL1THTTree = fs -> make<TTree>("genHTL1THTTree", "TTree with generator-level / L1T HT information");
    this -> _genHTL1THTTree -> Branch("genHT", &(this -> _genHT), "genHT/F");
    this -> _genHTL1THTTree -> Branch("l1tHT", &(this -> _l1tHT), "l1tHT/F");
  }
  //Creating a TTree where to save gen-l1t MHT, if I have received the MHT tag and the gen-jet tag required to compute it
  if ((this -> _genJetCollectionTag) && (this -> _l1tMHTCollectionTag))
  {
    this -> _genMHTL1TMHTTree = fs -> make<TTree>("genMHTL1TMHTTree", "TTree with generator-level / L1T MHT information");
    this -> _genMHTL1TMHTTree -> Branch("genMHT", &(this -> _genMHT), "genMHT/F");
    this -> _genMHTL1TMHTTree -> Branch("l1tMHT", &(this -> _l1tMHT), "l1tMHT/F");
  }
 
   
}


// obtains tokens to access data in the input CMSSW file if they are present

//the "try" statement tries to access the inputtag, if that fails the block in the catch statement gets executed
//I use this mechanism to create optional parameters

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
    this -> _l1tMHTCollectionTag = new edm::EDGetTokenT< BXVector<l1t::EtSum> >(consumes< BXVector<l1t::EtSum> > (iConfig.getParameter< edm::InputTag >("l1tMHTCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> MHT configuration not found, proceeding without adding MHT into to tree" << std::endl;
    this -> _l1tMHTCollectionTag = NULL;
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

//liberating the memory area held by the tokens
void SaveGenSumsAndL1Sums::_freeTokens()
{
  if (this -> _genMETCollectionTag) delete this -> _genMETCollectionTag;
  if (this -> _l1tMETCollectionTag) delete this -> _l1tMETCollectionTag;
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _l1tHTCollectionTag) delete this -> _l1tHTCollectionTag;
  if (this -> _l1tMHTCollectionTag) delete this -> _l1tMHTCollectionTag;
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
  //for each event I check if I am supposed to save MET by looking up if the tokens are properly initialised (i.e. not NULL)
  //if they are properly initialised I store the gen met and l1t met in the variables that will be than pushed to the TTree
  if ((this -> _genMETCollectionTag) && (this -> _l1tMETCollectionTag))
  {
    //token are initialised, I use them to access the data
    edm::Handle < std::vector< reco::GenMET> > lGenMETCollectionHandle;
    edm::Handle < BXVector<l1t::EtSum> > lL1TMETCollectionHandle;
    iEvent.getByToken(*(this -> _genMETCollectionTag), lGenMETCollectionHandle);
    iEvent.getByToken(*(this -> _l1tMETCollectionTag), lL1TMETCollectionHandle);
    // retrieving gen MET and saving it to the genMet memory area
    this -> _genMET = lGenMETCollectionHandle -> front().pt();
    // retrieving sums
    const BXVector<l1t::EtSum> & lL1TMETCollection = *lL1TMETCollectionHandle;
    // looking up the sum containing the MET and storing it to the l1tMet memory area
    for (const l1t::EtSum & lL1TMET: lL1TMETCollection) if (lL1TMET.getType() == l1t::EtSum::EtSumType::kMissingEt) this -> _l1tMET = lL1TMET.pt();
    // pushing the gen-l1t MET to the tree
    this -> _genMETL1TMETTree -> Fill();
  }

  //checking if the HT tokens are initialised (I require the l1tHT token and the genJet collection that I will use to compute the genHT)
  if ((this -> _genJetCollectionTag) && (this -> _l1tHTCollectionTag))
  {
    //tokens are initialised, accessing the data
    edm::Handle < std::vector< reco::GenJet> > lGenJetCollectionHandle;
    edm::Handle < BXVector<l1t::EtSum> > lL1THTCollectionHandle;
    iEvent.getByToken(*(this -> _genJetCollectionTag), lGenJetCollectionHandle);
    iEvent.getByToken(*(this -> _l1tHTCollectionTag), lL1THTCollectionHandle);
    // access l1t HT data 
    const BXVector<l1t::EtSum> & lL1THTCollection = *lL1THTCollectionHandle;
    // looking up the sum containing the l1tHET and storing it in the l1tMHT memory area
    for (const l1t::EtSum & lL1THT: lL1THTCollection) if (lL1THT.getType() == l1t::EtSum::EtSumType::kTotalHt) this -> _l1tHT = lL1THT.pt();
    // computing gen HT
    this -> _genHT = this -> _computeHT(*lGenJetCollectionHandle);
    //pushing gen-l1t ht to the tree
    this -> _genHTL1THTTree -> Fill();
  }

  //checking if the MHT tokens are initialised
  if ((this -> _genJetCollectionTag) && (this -> _l1tMHTCollectionTag))
  {
    edm::Handle < std::vector< reco::GenJet> > lGenJetCollectionHandle;
    edm::Handle < BXVector<l1t::EtSum> > lL1TMHTCollectionHandle;
    iEvent.getByToken(*(this -> _genJetCollectionTag), lGenJetCollectionHandle);
    iEvent.getByToken(*(this -> _l1tMHTCollectionTag), lL1TMHTCollectionHandle);
    // l1t MHT data
    const BXVector<l1t::EtSum> & lL1TMHTCollection = *lL1TMHTCollectionHandle;
    // storing in l1tMHT memory area
    for (const l1t::EtSum & lL1TMHT: lL1TMHTCollection) if (lL1TMHT.getType() == l1t::EtSum::EtSumType::kMissingHt) this -> _l1tMHT = lL1TMHT.pt();
    //computing gen MHT
    this -> _genMHT = this -> _computeMHT(*lGenJetCollectionHandle);

    // if((this -> _genMHT > 275)) std::cout << "l1tMHT: " << this -> _l1tMHT << "\tgenMHT: " << this -> _genMHT << std::endl;
    //pushing gen-l1t ht to the tree
    this -> _genMHTL1TMHTTree -> Fill();
  }

}

// computes the genHT starting from the gen jet collection
// i loop over every gen jet and sum the pt to the total HT of the event if their momentum is above the htPtThreshold
double SaveGenSumsAndL1Sums::_computeHT(const std::vector<reco::GenJet>& jetVector) 
{
  double lHT = 0;

  for (const auto & jet: jetVector)
  {
    lHT += (jet.pt() >= this -> _htPtThreshold) ? jet.pt() : 0;
  }

  return lHT;
}

// computes the genMHT starting from the gen jet collection
double SaveGenSumsAndL1Sums::_computeMHT(const std::vector<reco::GenJet>& jetVector) 
{
  double lTotalJetPx = 0;
  double lTotalJetPy = 0;
  
  for (const auto & jet: jetVector)
  {
    // checking if above threshold
    lTotalJetPx += (jet.pt() >= this -> _htPtThreshold) ? jet.px() : 0;
    lTotalJetPy += (jet.pt() >= this -> _htPtThreshold) ? jet.py() : 0;
  }
  double lMHT = sqrt(lTotalJetPx * lTotalJetPx + lTotalJetPy * lTotalJetPy);

  return lMHT;
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