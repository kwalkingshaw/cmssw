// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: ConvertGenJetToL1TJet
//
/**\class ConvertGenJetToL1TJet ConvertGenJetToL1TJet.cc L1Trigger/L1CaloTrigger/plugin/ConvertGenJetToL1TJet.cc

Description: Produces jets with sliding window algorithm using pfcluster and pfcandidates

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/L1Trigger/interface/Jet.h"


#include "TH2F.h"

#include <csignal>

//UNCOMMENT TO CREATE DEBUG HISTO
//#define DEBUG

//class ConvertGenJetToL1TJet : public edm::EDProducer {
class ConvertGenJetToL1TJet : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit ConvertGenJetToL1TJet(const edm::ParameterSet&);
      ~ConvertGenJetToL1TJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      edm::EDGetTokenT<std::vector<reco::GenJet>> *_ak4GenJetFromPfCandidatesCollectionTag;
      edm::EDGetTokenT<std::vector<reco::GenJet>> *_ak4GenJetFromPfClustersCollectionTag;
      BXVector<l1t::Jet> _convertGenJetCollectionToL1TCollection(const std::vector<reco::GenJet> &);

};

ConvertGenJetToL1TJet::ConvertGenJetToL1TJet(const edm::ParameterSet& iConfig)
{

  try
  {
    this -> _ak4GenJetFromPfCandidatesCollectionTag = new edm::EDGetTokenT< std::vector<reco::GenJet> >(consumes< std::vector<reco::GenJet> > (iConfig.getParameter< edm::InputTag >("ak4GenJetFromPfCandidatesCollectionTag")));
    produces<BXVector<l1t::Jet> >( "ak4L1TJetFromPfCandidates" ).setBranchAlias("ak4L1TJetFromPfCandidates");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> ak4GenJetFromPfCandidatesCollectionTag not found" << std::endl;
    this -> _ak4GenJetFromPfCandidatesCollectionTag = NULL;
  }

  try
  {
    this -> _ak4GenJetFromPfClustersCollectionTag = new edm::EDGetTokenT< std::vector<reco::GenJet> >(consumes< std::vector<reco::GenJet> > (iConfig.getParameter< edm::InputTag >("ak4GenJetFromPfClustersCollectionTag")));
    produces<BXVector<l1t::Jet> >( "ak4L1TJetFromPfClusters" ).setBranchAlias("ak4L1TJetFromPfClusters");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> ak4GenJetFromPfClustersCollectionTag not found" << std::endl;
    this -> _ak4GenJetFromPfClustersCollectionTag = NULL;
  }
  
}

ConvertGenJetToL1TJet::~ConvertGenJetToL1TJet()
{
  if (this -> _ak4GenJetFromPfCandidatesCollectionTag ) delete this -> _ak4GenJetFromPfCandidatesCollectionTag;
  if (this -> _ak4GenJetFromPfClustersCollectionTag ) delete this -> _ak4GenJetFromPfClustersCollectionTag;
}

void ConvertGenJetToL1TJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if (this -> _ak4GenJetFromPfCandidatesCollectionTag) {
    edm::Handle < std::vector< reco::GenJet > > ak4GenJetFromPfCandidatesCollectionHandle;
    iEvent.getByToken(*(this -> _ak4GenJetFromPfCandidatesCollectionTag), ak4GenJetFromPfCandidatesCollectionHandle);
    const auto l1jetVector = this -> _convertGenJetCollectionToL1TCollection(*ak4GenJetFromPfCandidatesCollectionHandle);
    std::unique_ptr< BXVector<l1t::Jet> > l1jetVectorPtr(new BXVector<l1t::Jet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), "ak4L1TJetFromPfCandidates");
  }

  if (this -> _ak4GenJetFromPfClustersCollectionTag) {
    edm::Handle < std::vector< reco::GenJet > > ak4GenJetFromPfClustersCollectionHandle;
    iEvent.getByToken(*(this -> _ak4GenJetFromPfClustersCollectionTag), ak4GenJetFromPfClustersCollectionHandle);
    const auto l1jetVector = this -> _convertGenJetCollectionToL1TCollection(*ak4GenJetFromPfClustersCollectionHandle);
    std::unique_ptr< BXVector<l1t::Jet> > l1jetVectorPtr(new BXVector<l1t::Jet>(l1jetVector));
    iEvent.put(std::move(l1jetVectorPtr), "ak4L1TJetFromPfClusters");
  }

  return;

}

BXVector<l1t::Jet> ConvertGenJetToL1TJet::_convertGenJetCollectionToL1TCollection(const std::vector<reco::GenJet> & genJets)
{

  BXVector<l1t::Jet> jets;
  for (const auto& genJet: genJets){
    math::PtEtaPhiMLorentzVector ptVector;
    ptVector.SetPt(genJet.pt());
    ptVector.SetEta(genJet.eta());
    ptVector.SetPhi(genJet.phi());

    l1t::Jet jet(ptVector);
    jets.push_back(0, jet);
  }

  return jets;
}

void
ConvertGenJetToL1TJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConvertGenJetToL1TJet);