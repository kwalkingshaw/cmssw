// -*- C++ -*-
//
// Package:    TriggerPerformanceAnalysis/StoreCandidatesToTree
// Class:      StoreCandidatesToTree
//
/**\class StoreCandidatesToTree StoreCandidatesToTree.cc TriggerPerformanceAnalysis/DataValidation/plugins/StoreCandidatesToTree.cc

 Description: Saves a list of candidates to root tree

 Implementation:
    Creates a branch in a TTree with an array of candidates. For each candidate pt, eta, phi is saved.
*/
//
// Original Author:  Simone Bologna
//         Created:  Tue, 18 Feb 2020 09:50:46 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TTree.h"


// class declaration

class StoreCandidatesToTree : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit StoreCandidatesToTree(const edm::ParameterSet&);
    ~StoreCandidatesToTree();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    edm::EDGetTokenT< edm::View< reco::Candidate> > _candidateCollectionTag;
    std::string _treeName;
    TTree* _candidateTree;
    unsigned int _maxNumberOfCandidates;
    unsigned int _numberOfCandidates;
    float* _candidatePt;
    float* _candidateEta;
    float* _candidatePhi;
    // ----------member data ---------------------------
};

//
// constructors and destructor
//
// retrieves the tags for the candidate list
// prepares the data holders
// creates a tree and branches where the candidate list is going to be stored
StoreCandidatesToTree::StoreCandidatesToTree(const edm::ParameterSet& iConfig):
  _candidateCollectionTag(consumes< edm::View< reco::Candidate > >  (iConfig.getParameter< edm::InputTag>  ("candidateCollectionTag"))),
  _treeName(iConfig.getParameter< std::string> ("treeName")),
  _maxNumberOfCandidates(iConfig.getParameter< unsigned int> ("maxNumberOfCandidates"))
{
  
  //preparing data holders
  this -> _candidatePt = new float[this -> _maxNumberOfCandidates];
  this -> _candidateEta = new float[this -> _maxNumberOfCandidates];
  this -> _candidatePhi = new float[this -> _maxNumberOfCandidates];

  // Accessing the output file via a TFileService
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  // creating the tree in the file
  this -> _candidateTree = fs -> make<TTree>(this -> _treeName.c_str(), "Tree with list of candidates");
  // adding the branches
  this -> _candidateTree -> Branch("length", &(this -> _numberOfCandidates), "length/i");
  this -> _candidateTree -> Branch("pt", this -> _candidatePt, "pt[length]/F");
  this -> _candidateTree -> Branch("eta", this -> _candidateEta, "eta[length]/F");
  this -> _candidateTree -> Branch("phi", this -> _candidatePhi, "phi[length]/F");

}

// deallocating memory used by data holder
StoreCandidatesToTree::~StoreCandidatesToTree()
{
  delete this -> _candidatePt;
  delete this -> _candidateEta;
  delete this -> _candidatePhi;
}


//
// member functions
//

// ------------ method called for each event  ------------
// accesses candidates, saves data to data holders, fills the tree
void
StoreCandidatesToTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //accessing candidates
  edm::Handle < edm::View< reco::Candidate > > candidateCollectionHandle;
  iEvent.getByToken(this -> _candidateCollectionTag, candidateCollectionHandle);

  //filling data holders
  this -> _numberOfCandidates = candidateCollectionHandle -> size();
  for (unsigned int x = 0; x < this -> _numberOfCandidates; x++)
  {
    const reco::Candidate & lCandidate = candidateCollectionHandle -> at(x);
    this -> _candidatePt[x] =  lCandidate.pt();
    this -> _candidateEta[x] =  lCandidate.eta();
    this -> _candidatePhi[x] =  lCandidate.phi();
    
    // dumps candidates for debuggin'
    // std::cout << this -> _candidatePt[x] << "\t" << this -> _candidateEta[x] << "\t" << this -> _candidatePhi[x] << std::endl;
  }

  this -> _candidateTree -> Fill();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
StoreCandidatesToTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(StoreCandidatesToTree);
