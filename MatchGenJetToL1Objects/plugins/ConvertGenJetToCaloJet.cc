// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: ConvertGenJetToCaloJet
//
/**\class ConvertGenJetToCaloJet ConvertGenJetToCaloJet.cc L1Trigger/L1CaloTrigger/plugin/ConvertGenJetToCaloJet.cc

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
#include "DataFormats/JetReco/interface/CaloJet.h"
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

//class ConvertGenJetToCaloJet : public edm::EDProducer {
class ConvertGenJetToCaloJet : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit ConvertGenJetToCaloJet(const edm::ParameterSet&);
      ~ConvertGenJetToCaloJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      edm::EDGetTokenT<std::vector<reco::GenJet>> *_ak4GenJetFromPfCandidatesCollectionTag;
      edm::EDGetTokenT<std::vector<reco::GenJet>> *_ak4GenJetFromPfClustersCollectionTag;
      std::vector<reco::CaloJet> _convertGenJetCollectionToCaloJetCollection(const std::vector<reco::GenJet> &);

};

ConvertGenJetToCaloJet::ConvertGenJetToCaloJet(const edm::ParameterSet& iConfig)
{

  try
  {
    this -> _ak4GenJetFromPfCandidatesCollectionTag = new edm::EDGetTokenT< std::vector<reco::GenJet> >(consumes< std::vector<reco::GenJet> > (iConfig.getParameter< edm::InputTag >("ak4GenJetFromPfCandidatesCollectionTag")));
    produces<std::vector<reco::CaloJet> >( "ak4CaloJetFromPfCandidates" ).setBranchAlias("ak4CaloJetFromPfCandidates");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> ak4GenJetFromPfCandidatesCollectionTag not found" << std::endl;
    this -> _ak4GenJetFromPfCandidatesCollectionTag = NULL;
  }

  try
  {
    this -> _ak4GenJetFromPfClustersCollectionTag = new edm::EDGetTokenT< std::vector<reco::GenJet> >(consumes< std::vector<reco::GenJet> > (iConfig.getParameter< edm::InputTag >("ak4GenJetFromPfClustersCollectionTag")));
    produces<std::vector<reco::CaloJet> >( "ak4CaloJetFromPfClusters" ).setBranchAlias("ak4CaloJetFromPfClusters");
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> ak4GenJetFromPfClustersCollectionTag not found" << std::endl;
    this -> _ak4GenJetFromPfClustersCollectionTag = NULL;
  }
  
}

ConvertGenJetToCaloJet::~ConvertGenJetToCaloJet()
{
  if (this -> _ak4GenJetFromPfCandidatesCollectionTag ) delete this -> _ak4GenJetFromPfCandidatesCollectionTag;
  if (this -> _ak4GenJetFromPfClustersCollectionTag ) delete this -> _ak4GenJetFromPfClustersCollectionTag;
}

void ConvertGenJetToCaloJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if (this -> _ak4GenJetFromPfCandidatesCollectionTag) {
    edm::Handle < std::vector< reco::GenJet > > ak4GenJetFromPfCandidatesCollectionHandle;
    iEvent.getByToken(*(this -> _ak4GenJetFromPfCandidatesCollectionTag), ak4GenJetFromPfCandidatesCollectionHandle);
    const auto calojetVector = this -> _convertGenJetCollectionToCaloJetCollection(*ak4GenJetFromPfCandidatesCollectionHandle);
    std::unique_ptr< std::vector<reco::CaloJet> > calojetVectorPtr(new std::vector<reco::CaloJet>(calojetVector));
    iEvent.put(std::move(calojetVectorPtr), "ak4CaloJetFromPfCandidates");
  }

  if (this -> _ak4GenJetFromPfClustersCollectionTag) {
    edm::Handle < std::vector< reco::GenJet > > ak4GenJetFromPfClustersCollectionHandle;
    iEvent.getByToken(*(this -> _ak4GenJetFromPfClustersCollectionTag), ak4GenJetFromPfClustersCollectionHandle);
    const auto calojetVector = this -> _convertGenJetCollectionToCaloJetCollection(*ak4GenJetFromPfClustersCollectionHandle);
    std::unique_ptr< std::vector<reco::CaloJet> > calojetVectorPtr(new std::vector<reco::CaloJet>(calojetVector));
    iEvent.put(std::move(calojetVectorPtr), "ak4CaloJetFromPfClusters");
  }

  return;

}

std::vector<reco::CaloJet> ConvertGenJetToCaloJet::_convertGenJetCollectionToCaloJetCollection(const std::vector<reco::GenJet> & genJets)
{

  std::vector<reco::CaloJet> jets;
  for (const auto& genJet: genJets){
    math::PtEtaPhiMLorentzVector ptVector;
    ptVector.SetPt(genJet.pt());
    ptVector.SetEta(genJet.eta());
    ptVector.SetPhi(genJet.phi());

    reco::CaloJet jet;
    jet.setP4(ptVector);
    jets.push_back(jet);
  }

  return jets;
}

void
ConvertGenJetToCaloJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConvertGenJetToCaloJet);