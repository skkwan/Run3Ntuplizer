#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "L1Trigger/Run3Ntuplizer/interface/Run3Ntuplizer.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>

using namespace edm;
using std::cout;
using std::endl;
using std::vector;

bool compareByPtJets (l1extra::L1JetParticle i,l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };
bool compareByPtTaus (l1t::Tau i,l1t::Tau j) { return(i.pt()>j.pt()); };

Run3Ntuplizer::Run3Ntuplizer( const ParameterSet & cfg ) :
  genSrc_((cfg.getParameter<edm::InputTag>( "genParticles"))),
  tauSrc_(       consumes<vector<pat::Tau>     >(cfg.getParameter<edm::InputTag>("miniTaus"))),
  stage2TauSrc_( consumes<vector <l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("stage2Taus" ))),
  stage2IsoTauSrc_( consumes<vector<l1extra::L1JetParticle> >(cfg.getParameter<edm::InputTag>("stage2IsoTaus"))),
  stage2DigisTauSrc_( consumes<BXVector<l1t::Tau> >(cfg.getParameter<edm::InputTag>("stage2DigisTaus")))
{
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  
  folderName_          = cfg.getUntrackedParameter<std::string>("folderName");
  isData_              = cfg.getParameter<bool>("isData");
  folder               = tfs_->mkdir(folderName_);
  
  efficiencyTree = folder.make<TTree>("efficiencyTree", "Gen Matched Tau Tree");
  createBranchesTau(efficiencyTree);

}

void Run3Ntuplizer::createBranchesTau(TTree *tree){

  tree->Branch("l1TauPt",       &l1TauPt,   "l1TauPt/D");
  tree->Branch("l1TauEta",      &l1TauEta,  "l1TauEta/D");
  tree->Branch("l1TauPhi",      &l1TauPhi,  "l1TauPhi/D");

  tree->Branch("recoTauPt",     &recoTauPt, "recoTauPt/D");
  tree->Branch("recoTauEta",    &recoTauEta,"recoTauEta/D");
  tree->Branch("recoTauPhi",    &recoTauPhi,"recoTauPhi/D");
  tree->Branch("recoTauDM",     &recoTauDM, "recoTauDM/D");

  tree->Branch("genTauPt",      &genTauPt,  "genTauPt/D");
  tree->Branch("genTauEta",     &genTauEta, "genTauEta/D");
  tree->Branch("genTauPhi",     &genTauPhi, "genTauPhi/D");
  tree->Branch("genTauDM",      &genTauDM,  "genTauDM/D");


}

void Run3Ntuplizer::beginJob( const EventSetup & es) {
  cout<<"begin job..."<<std::endl;
}

void Run3Ntuplizer::analyze( const Event& evt, const EventSetup& es )
{
  cout<<"Analyzing..."<<std::endl;

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  edm::Handle < vector<l1extra::L1JetParticle> > stage2Taus;
  edm::Handle < vector<l1extra::L1JetParticle> > stage2IsoTaus;
  edm::Handle < BXVector<l1t::Tau> > stage2DigiTaus;

  edm::Handle< std::vector<pat::Tau> > miniTaus;

  if(!evt.getByToken( tauSrc_, miniTaus))
    cout<<"No miniAOD particles found"<<std::endl;

  if(!evt.getByToken(stage2TauSrc_, stage2Taus))
    cout<<"ERROR GETTING THE STAGE 2 TAUS"<<std::endl;
  else
    cout<<"Stage2 Tau Size: "<<stage2Taus->size()<<std::endl;

  if(!evt.getByToken(stage2DigisTauSrc_, stage2DigiTaus))
    cout<<"ERROR GETTING THE STAGE 2 TAUS"<<std::endl;
  else
    cout<<"Stage2 Digi Tau Size: "<<stage2DigiTaus->size()<<std::endl;

  if(!evt.getByToken(stage2IsoTauSrc_, stage2IsoTaus))
    cout<<"ERROR GETTING THE STAGE 2 ISO TAUS"<<std::endl;

  //sort the L1 taus
  vector<l1t::Tau> l1TausSorted;
  vector<l1extra::L1JetParticle> l1IsoTausSorted;
  for( BXVector<l1t::Tau>::const_iterator l1Tau = stage2DigiTaus->begin(); l1Tau != stage2DigiTaus->end(); l1Tau++ ){
    l1TausSorted.push_back(*l1Tau);
  }

  for( vector<l1extra::L1JetParticle>::const_iterator l1Tau = stage2IsoTaus->begin(); l1Tau != stage2IsoTaus->end(); l1Tau++ ){
    l1IsoTausSorted.push_back(*l1Tau);
    if(abs(l1Tau->eta()) < 2.4) l1IsoTausSorted.push_back(*l1Tau);
    cout<<"l1Tau Pt: "<<l1Tau->pt()<<" Eta: "<<l1Tau->eta()<<" Phi: "<<l1Tau->phi()<<std::endl;
  }

  std::sort(l1TausSorted.begin(),l1TausSorted.end(),compareByPtTaus);
  std::sort(l1IsoTausSorted.begin(),l1IsoTausSorted.end(),compareByPtJets);

  // Now for the Taus
  edm::Handle<GenParticleCollectionType> genParticleHandle;
  if(!isData_){
    if(!evt.getByToken(genToken_,genParticleHandle))
      cout<<"No gen Particles Found "<<std::endl;
  }  

  vector<reco::GenParticle> genTaus;
  vector<reco::GenParticle> genParticles;
  vector<genVisTau> genVisTaus;
  genVisTaus.clear();

  if(!isData_){
    for(unsigned int i = 0; i< genParticleHandle->size(); i++){
      edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
      genParticles.push_back(*ptr);
      /*
      if(abs(ptr->pdgId())==111 && abs(ptr->eta()<1.74)){
      genPiZeros.push_back(*ptr);
      //cout<<"Found PiZero PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
      }
      if(abs(ptr->pdgId())==211 && abs(ptr->eta()<1.74)){
      genPiPluss.push_back(*ptr);
      //cout<<"Found PiPlus PDGID 111 pt: "<<ptr->pt()<<" eta: "<<ptr->eta()<<" phi: "<<ptr->phi()<<std::endl;
      }*/
      if(abs(ptr->pdgId())==15){
	genTaus.push_back(*ptr);
      }
    }

    for(auto genTau: genTaus){
      reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &genParticles);
      genVisTau Temp;
      int decayMode = GetDecayMode(&genTau);
      Temp.p4 = visGenTau;
      Temp.decayMode = decayMode;
      genVisTaus.push_back(Temp);
      //cout<<"Tau Decay Mode "<<decayMode<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;
    }
  }

  //  zeroOutAllVariables();
  
  if(miniTaus->size() > 0){
    recoTauPt        = -99;
    recoTauEta       = -99;
    recoTauPhi       = -99;
    //recoChargedIso = -99;
    //recoNeutralIso = -99;
    //recoRawIso     = -99;
    recoTauDM  = -99;
    
    if(miniTaus->at(0).tauID("decayModeFinding")>0){
      recoTauPt         = miniTaus->at(0).p4().Pt();
      recoTauEta        = miniTaus->at(0).p4().Eta();
      recoTauPhi        = miniTaus->at(0).p4().Phi();
      //recoChargedIso = miniTaus->at(i).tauID("chargedIsoPtSum");
      //recoNeutralIso = miniTaus->at(i).tauID("neutralIsoPtSum");
      //recoRawIso     = miniTaus->at(i).tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      recoTauDM  = miniTaus->at(0).decayMode();
      //cout<<"=== Found recoTau: "<<recoPt<<" Eta: "<<recoEta<<" Phi: "<< recoPhi <<std::endl;
      
      l1TauPt  = -99;
      l1TauEta = -99;
      l1TauPhi = -99;
      
      for(unsigned int i = 0; i < l1TausSorted.size(); i++){
	if(( reco::deltaR(l1TausSorted.at(i).eta(), l1TausSorted.at(i).phi(), 
			  recoTauEta, recoTauPhi) < 0.5 )
	   && (l1TausSorted.at(i).pt() > l1TauPt))
	  {
	    l1TauPt  = l1TausSorted.at(i).pt();
	    l1TauEta = l1TausSorted.at(i).eta();
	    l1TauPhi = l1TausSorted.at(i).phi();
	    break;
	  }
      }

      genTauPt  = -99;
      genTauEta = -99;
      genTauPhi = -99;
      genTauDM  = -99;

      if(!isData_)
	for(auto genTau :genVisTaus){
	  if(( reco::deltaR(genTau.p4.eta(), genTau.p4.phi(), 
			    recoTauEta, recoTauPhi) < 0.5 )){
	    genTauPt  = genTau.p4.pt();
	    genTauEta = genTau.p4.eta();
	    genTauPhi = genTau.p4.phi();
	    genTauDM  = genTau.decayMode;
	    break;
	  }
	}
      
    }
  }
  
  efficiencyTree->Fill();

}


void Run3Ntuplizer::zeroOutAllVariables(){

  recoTauPt = -99;  recoTauEta = -99;  recoTauPhi = -99;  
  l1TauPt   = -99;  l1TauEta   = -99;  recoTauEta = -99;


};

void Run3Ntuplizer::endJob() {
}

Run3Ntuplizer::~Run3Ntuplizer(){

}

DEFINE_FWK_MODULE(Run3Ntuplizer);
