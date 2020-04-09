#ifndef Run3Ntuplizer_H
#define Run3Ntuplizer_H

// system include files
#include <memory>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"
// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "L1Trigger/L1TCaloLayer1/src/L1UCTCollections.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"

/* TMVA */
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

//
// class declaration
//
using std::vector;

class Run3Ntuplizer : public edm::EDAnalyzer {

 public:
  
  // Constructor
  Run3Ntuplizer(const edm::ParameterSet& ps);
  
  // Destructor
  virtual ~Run3Ntuplizer();

  edm::Service<TFileService> tfs_;

  // File I/O
  std::ofstream file0, file1, file10;
  TFileDirectory folder;
  
  // Trees
  TTree* efficiencyTree;

  // Other variables
  int run, lumi, event;
  float nvtx;

  // Helper functions
  void zeroOutAllVariables();
  void createBranchesTau(TTree *tree);

 protected:
  // Analyze
  void analyze(const edm::Event& evt, const edm::EventSetup& es);
  
  // BeginJob
  void beginJob(const edm::EventSetup &es);
  
  // EndJob
  void endJob(void);

 private:
  // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;
  int nev_; // Number of events processed
  bool verbose_;
  std::ofstream logFile_;
  edm::InputTag rctSource_; 
  edm::InputTag genSrc_;

  edm::EDGetTokenT<vector<pat::Tau> > tauSrc_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > stage2TauSrc_;
  edm::EDGetTokenT<vector <l1extra::L1JetParticle> > stage2IsoTauSrc_;
  edm::EDGetTokenT<BXVector<l1t::Tau> > stage2DigisTauSrc_;
  
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;

  std::string folderName_;

  double genTauPt,  genTauEta,  genTauPhi,  genTauDM;
  double recoTauPt, recoTauEta, recoTauPhi, recoTauDM;
  double recoRawIso, recoChargedIso, recoNeutralIso;
  double l1TauPt, l1TauEta, l1TauPhi;
  int l1IsoEt;

  bool isData_;

  

  
};

#endif
