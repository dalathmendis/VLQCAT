// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      OS2LAna
// 
/**\class VLQAna OS2LAna.cc Analysis/VLQAna/plugins/OS2LAna.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
// Modified: Sadia Khalil
//           25 Mar 2016 17:11 CDT
//

#include <iostream>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateFilter.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ApplyLeptonIDSFs.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/METMaker.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"
#include "Analysis/VLQAna/interface/ApplyLeptonTrigSFs.h"
#include "Analysis/VLQAna/interface/DYNLOEwkKfact.h"
#include "Analysis/VLQAna/interface/OS2LTree.h"

#include "Analysis/VLQAna/interface/TopCandsProducer.h"
#include "Analysis/VLQAna/interface/ZCandsProducer.h"
#include "Analysis/VLQAna/interface/HCandsProducer.h"
#include "Analysis/VLQAna/interface/ZHCandsProducer.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

//
// class declaration
//

class OS2LAna : public edm::EDFilter {
  public:
    explicit OS2LAna(const edm::ParameterSet&);
    ~OS2LAna();

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt);
    // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<double>   t_evtwtPVLow      ;
    edm::EDGetTokenT<double>   t_evtwtPVHigh     ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision     ;
    edm::EDGetTokenT<int>      t_evtno           ; 
    edm::EDGetTokenT<int>      t_runno           ;
    edm::EDGetTokenT<int>      t_lumisec         ;
    edm::ParameterSet DilepCandParams_           ; 
    edm::ParameterSet ZCandParams_               ; 
    //edm::ParameterSet BoostedZCandParams_        ; 
    //edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const unsigned int NAK4Min_                  ;
    const double HTMin_                          ;
    const double STMin_                          ; 
    const double STMaxControl_                   ;
    const bool skim_                             ;
    const bool filterSignal_                     ;
    const bool additionalPlots_                  ;
    const std::string signalType_                ;
    const std::string zdecayMode_                ;
    const bool btagsf_bcUp_                      ;
    const bool btagsf_bcDown_                    ;
    const bool btagsf_lUp_                       ;
    const bool btagsf_lDown_                     ;
    const bool PileupUp_                         ;
    const bool PileupDown_                       ;   
    const bool categorize_                      ;  
    const bool applyLeptonIDSFs_                 ;
    const bool applyLeptonTrigSFs_               ;
    const bool applyBTagSFs_                     ;
    const bool applyDYNLOCorr_                   ;
    const std::string fname_DYNLOCorr_           ; 
    const std::string funname_DYNLOCorr_         ; 
    DYNLOEwkKfact dynloewkkfact                  ;
    ApplyLeptonIDSFs lepIdSFs                    ;
    ApplyLeptonTrigSFs lepTrigSFs                ;
    METMaker metmaker                            ;
    MuonMaker muonmaker                          ; 
    ElectronMaker electronmaker                  ; 
    JetMaker jetAK4maker                         ; 
    JetMaker jetAK4BTaggedmaker                  ; 
    JetMaker jetAK8maker                         ; 
    JetMaker jetHTaggedmaker                     ; 
    JetMaker jetWTaggedmaker                     ; 
    JetMaker jetTopTaggedmaker                   ; 
    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 
    PickGenPart genpart                          ;
    const std::string fnamebtagSF_               ;
    const std::string fnameSJbtagSF_             ;
    const std::string btageffmap_                ;
    std::unique_ptr<BTagSFUtils> btagsfutils_    ; 
    std::unique_ptr<BTagSFUtils> sjbtagsfutils_  ; 

    const bool maketree_ ;
    TTree* tree_ ; 
    os2l::OS2LAnaTree os2ltree_ ; 
};

using namespace std;

// static data member definitions
void OS2LAna::fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt){

  for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
    float scEta = goodElectrons.at(iele).getscEta();
    if(fabs(scEta) <= 1.479){
      h1_["Eta_EB_el_pre"]-> Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EB_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EB_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EB_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EB_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EB_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
      h1_["SCETA_EB_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EB_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EB_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EB_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EB_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EB_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);


    }
    else if  (fabs(scEta) > 1.479 && fabs(scEta) < 2.4){
      h1_["Eta_EE_el_pre"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EE_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EE_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EE_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EE_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
      h1_["SCETA_EE_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EE_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EE_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EE_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EE_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EE_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);

    }
  }
}

// constructors and destructor
OS2LAna::OS2LAna(const edm::ParameterSet& iConfig) :
  t_evttype               (consumes<string>  (iConfig.getParameter<edm::InputTag>("evttype"))),
  t_evtwtGen              (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
  t_evtwtPV               (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
  t_evtwtPVLow            (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVLow"))),
  t_evtwtPVHigh           (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),
  t_npv                   (consumes<unsigned>(iConfig.getParameter<edm::InputTag>("npv"))),
  t_hltdecision           (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecision"))),
  t_evtno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("evtno"))),
  t_runno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("runno"))),
  t_lumisec               (consumes<int>     (iConfig.getParameter<edm::InputTag>("lumisec"))),
  DilepCandParams_        (iConfig.getParameter<edm::ParameterSet> ("DilepCandParams")),
  ZCandParams_            (iConfig.getParameter<edm::ParameterSet> ("ZCandParams")),
  //BoostedZCandParams_     (iConfig.getParameter<edm::ParameterSet> ("BoostedZCandParams")),
  //GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  genParams_              (iConfig.getParameter<edm::ParameterSet> ("genParams")),
  NAK4Min_                (iConfig.getParameter<unsigned int>      ("NAK4Min")),
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  STMin_                  (iConfig.getParameter<double>            ("STMin")), 
  STMaxControl_           (iConfig.getParameter<double>            ("STMaxControl")), 
  skim_                   (iConfig.getParameter<bool>              ("skim")), 
  filterSignal_           (iConfig.getParameter<bool>              ("filterSignal")), 
  additionalPlots_        (iConfig.getParameter<bool>              ("additionalPlots")), 
  signalType_             (iConfig.getParameter<std::string>       ("signalType")), 
  zdecayMode_             (iConfig.getParameter<std::string>       ("zdecayMode")),
  btagsf_bcUp_            (iConfig.getParameter<bool>              ("btagsf_bcUp")),
  btagsf_bcDown_          (iConfig.getParameter<bool>              ("btagsf_bcDown")),
  btagsf_lUp_             (iConfig.getParameter<bool>              ("btagsf_lUp")),
  btagsf_lDown_           (iConfig.getParameter<bool>              ("btagsf_lDown")),
  PileupUp_               (iConfig.getParameter<bool>              ("PileupUp")),
  PileupDown_             (iConfig.getParameter<bool>              ("PileupDown")),
  categorize_             (iConfig.getParameter<bool>              ("categorize")),
  applyLeptonIDSFs_       (iConfig.getParameter<bool>              ("applyLeptonSFs")), 
  applyLeptonTrigSFs_     (iConfig.getParameter<bool>              ("applyLeptonTrigSFs")),
  applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")),
  applyDYNLOCorr_         (iConfig.getParameter<bool>              ("applyDYNLOCorr")),
  fname_DYNLOCorr_        (iConfig.getParameter<std::string>       ("File_DYNLOCorr")),
  funname_DYNLOCorr_      (iConfig.getParameter<std::string>       ("Fun_DYNLOCorr")),
  dynloewkkfact           (DYNLOEwkKfact(fname_DYNLOCorr_,funname_DYNLOCorr_)),
  lepIdSFs                (iConfig.getParameter<edm::ParameterSet> ("lepIdSFsParams")),
  lepTrigSFs              (iConfig.getParameter<edm::ParameterSet> ("lepTrigSFsParams")), 
  metmaker                (iConfig.getParameter<edm::ParameterSet> ("metselParams"),consumesCollector()),
  muonmaker               (iConfig.getParameter<edm::ParameterSet> ("muselParams"),consumesCollector()),
  electronmaker           (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector()),
  jetAK4maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4selParams"),consumesCollector()),
  jetAK4BTaggedmaker      (iConfig.getParameter<edm::ParameterSet> ("jetAK4BTaggedselParams"),consumesCollector()),
  jetAK8maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()),
  jetHTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetHTaggedselParams"),consumesCollector()),
  jetWTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetWTaggedselParams"),consumesCollector()),
  jetTopTaggedmaker       (iConfig.getParameter<edm::ParameterSet> ("jetTopTaggedselParams"),consumesCollector()),   
  genpart                 (genParams_, consumesCollector()),
  fnamebtagSF_            (iConfig.getParameter<std::string>       ("fnamebtagSF")),
  fnameSJbtagSF_          (iConfig.getParameter<std::string>       ("fnameSJbtagSF")),
  btageffmap_             (iConfig.getParameter<std::string>       ("btageffmap")),
  btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,20., 1000., 20., 1000., 20., 1000.,btageffmap_)),
  sjbtagsfutils_          (new BTagSFUtils(fnameSJbtagSF_,BTagEntry::OP_LOOSE,30., 450., 30., 450., 20., 1000.,btageffmap_)),
  maketree_               (iConfig.getParameter<bool>("maketree"))

{
  produces<vlq::JetCollection>("tjets") ; 
  produces<vlq::JetCollection>("wjets") ; 
  produces<vlq::JetCollection>("hjets") ;
  produces<vlq::JetCollection>("bjets") ; 
  produces<vlq::JetCollection>("jets") ; 
  produces<vlq::CandidateCollection>("zllcands") ; 
  produces<double>("PreWeight");
  produces<double>("btagsf");
  produces<double>("btagsfbcUp");
  produces<double>("btagsfbcDown");
  produces<double>("btagsflUp");
  produces<double>("btagsflDown");
  produces<double>("sjbtagsf");
  produces<double>("sjbtagsfbcUp");
  produces<double>("sjbtagsfbcDown");
  produces<double>("sjbtagsflUp");
  produces<double>("sjbtagsflDown");
  produces<double>("finalWeight");
}


OS2LAna::~OS2LAna() {}

bool OS2LAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<string>   h_evttype     ; evt.getByToken(t_evttype     ,h_evttype    ) ; 
  Handle<double>   h_evtwtGen    ; evt.getByToken(t_evtwtGen    ,h_evtwtGen   ) ; 
  Handle<double>   h_evtwtPV     ; evt.getByToken(t_evtwtPV     ,h_evtwtPV    ) ; 
  Handle<double>   h_evtwtPVLow  ; evt.getByToken(t_evtwtPVLow  ,h_evtwtPVLow ) ;
  Handle<double>   h_evtwtPVHigh ; evt.getByToken(t_evtwtPVHigh ,h_evtwtPVHigh) ;  
  Handle<unsigned> h_npv         ; evt.getByToken(t_npv         ,h_npv        ) ; 
  Handle<bool>     h_hltdecision ; evt.getByToken(t_hltdecision ,h_hltdecision) ; 
  Handle<int>      h_evtno       ; evt.getByToken(t_evtno       ,h_evtno      ) ; 
  Handle<int>      h_runno       ; evt.getByToken(t_runno       ,h_runno      ) ; 
  Handle<int>      h_lumisec     ; evt.getByToken(t_lumisec     ,h_lumisec    ) ; 

  const int npv(*h_npv.product());
  const int evtno(*h_evtno.product()) ;
  const int runno(*h_runno.product()) ;
  const int lumisec(*h_lumisec.product()) ;
  const bool isData(evtno > 0 ? true : false) ; 

  int signalType(-1);
  if (filterSignal_) {
    if( (skim_ || maketree_ ) && signalType_.empty() ){ 
      if      (*h_evttype.product() == "EvtType_MC_bZbZ") signalType = 1; 
      else if (*h_evttype.product() == "EvtType_MC_bZbH") signalType = 2; 
      else if (*h_evttype.product() == "EvtType_MC_bZtW") signalType = 3; 
      else if (*h_evttype.product() == "EvtType_MC_bHbH") signalType = 4; 
      else if (*h_evttype.product() == "EvtType_MC_bHtW") signalType = 5; 
      else if (*h_evttype.product() == "EvtType_MC_tWtW") signalType = 6; 
      else if (*h_evttype.product() == "EvtType_MC_tZtZ") signalType = 7; 
      else if (*h_evttype.product() == "EvtType_MC_tZtH") signalType = 8; 
      else if (*h_evttype.product() == "EvtType_MC_tZbW") signalType = 9; 
      else if (*h_evttype.product() == "EvtType_MC_tHtH") signalType = 10; 
      else if (*h_evttype.product() == "EvtType_MC_tHbW") signalType = 11; 
      else if (*h_evttype.product() == "EvtType_MC_bWbW") signalType = 12; 
      h1_["signalEvts_all"] -> Fill(signalType);
    }
    else{
      if (*h_evttype.product()!=signalType_) return false ;
      else  h1_["signalEvts"] -> Fill(1) ;
    }
  }

  double evtwt(1.0);
  if (PileupUp_)         evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVHigh.product()) ; 
  else if (PileupDown_)  evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVLow.product()) ;
  else                   evtwt = (*h_evtwtGen.product()) * (*h_evtwtPV.product()) ;

  h1_["cutflow"] -> Fill(1, evtwt) ;

  const bool hltdecision(*h_hltdecision.product()) ; 
  if ( hltdecision ) h1_["cutflow"] -> Fill(2, evtwt) ;
  else return false; //// Presel: HLT  

  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MetCollection goodMet;
  metmaker(evt, goodMet) ;
  vlq::CandidateCollection tops, W, B, BC, D, Z, ZB, H, Hb, ZH, ZHb;
  vlq::CandidateCollection dimuons, dielectrons, dileptons;   
  vlq::CandidateCollection zll; //generic collection

  // dilepton properties: M > 50, lead pt > 45, second pt > 25
  DileptonCandsProducer dileptonsprod(DilepCandParams_) ; 
  dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons); 
  dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons) ; 

  //================================================================
  //First pre-selection: 1) 2 OS dileptons from boosted Z, >=3 jets
  //================================================================

  //// Dilepton candidate
  if (zdecayMode_ == "zmumu") {dileptons = dimuons; }
  else if (zdecayMode_ == "zelel") {dileptons = dielectrons;}

  if (dileptons.size() > 0) h1_["cutflow"] -> Fill(3, evtwt) ;
  else return false ; //// Presel: dileptons

  //// Get Dy EWK correction
  if ( applyDYNLOCorr_ ) {
    double EWKNLOkfact( dynloewkkfact(dileptons.at(0).getPt()) ) ; ////GetDYNLOCorr(dileptons.at(0).getPt())) ; 
    evtwt *= EWKNLOkfact ;
  }

  //// Get lepton ID and Iso SF
  if  ( !isData ) {
    //    cout << "apply id sF " << applyLeptonIDSFs_ <<endl;
    if (applyLeptonIDSFs_) {
      if ( zdecayMode_ == "zmumu" ) {
        evtwt *= lepIdSFs.IDSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodMuons.at(1).getPt(),goodMuons.at(1).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(1).getPt(), goodMuons.at(1).getEta()) ; 
      }
      else if ( zdecayMode_ == "zelel" ) {
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(1).getPt(), goodElectrons.at(1).getEta()) ;
      }
    }

    if (applyLeptonTrigSFs_) {
      if ( zdecayMode_ == "zmumu" ) evtwt *= lepTrigSFs(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ; 
      else if ( zdecayMode_ == "zelel" ) evtwt *= lepTrigSFs(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ; 
    }

  } 

  //// Z mass candidate filter: 75 < M < 105, lead pt > 45, 2nd pt > 25, Z pt > 100
  CandidateFilter zllfilter(ZCandParams_) ; 
  zllfilter(dileptons, zll);

  //// jets
  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;

  // b-tagging:
  vlq::JetCollection goodBTaggedAK4Jets;
  jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ; 

  //// jet cleaning w.r.t dileptons
  CandidateCleaner cleanjets(0.4, -1); //// The second argument is for lepton 2D iso, setting to -1 disables it
  if (zdecayMode_ == "zmumu") {
    cleanjets(goodAK4Jets, goodMuons);
    cleanjets(goodBTaggedAK4Jets, goodMuons);
  }
  else if (zdecayMode_ == "zelel") {
    cleanjets(goodAK4Jets, goodElectrons);
    cleanjets(goodBTaggedAK4Jets, goodElectrons);
  } 

  //// At least one Z cand in event
  if(zll.size() == 1) {h1_["cutflow"] -> Fill(4, evtwt) ;}
  else return false ; //// Presel : Z->ll

  //// At least HT > 200 in event
  HT htak4(goodAK4Jets) ; 
  if ( htak4.getHT() > HTMin_ ) h1_["cutflow"] -> Fill(5, evtwt) ;  
  else return false ; //// Presel: HT cut

  //// At least 3 AK4 jets in event
  if (goodAK4Jets.size()  >= NAK4Min_ ) {h1_["cutflow"] -> Fill(6, evtwt) ;} 
  else return false; //// Presel N(AK4) 

  vlq::JetCollection goodAK8Jets ;
  jetAK8maker(evt, goodAK8Jets); 
  cleanjets(goodAK8Jets, goodMuons); 
  cleanjets(goodAK8Jets, goodElectrons); 

  double presel_wt(evtwt);
  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodAK4Jets) {
      csvs.push_back(jet.getCSV()) ; 
      pts.push_back(jet.getPt()) ; 
      etas.push_back(jet.getEta()) ; 
      flhads.push_back(jet.getHadronFlavour()) ; 
    }

    btagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetAK4BTaggedmaker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ;

    //// bTag SF along with sys. unc. options
    if (btagsf_bcUp_)
      evtwt *= btagsf_bcUp;
    else if (btagsf_bcDown_)
      evtwt *= btagsf_bcDown;
    else  if (btagsf_lUp_)
      evtwt *= btagsf_lUp;
    else if (btagsf_lDown_)
      evtwt *= btagsf_lDown;
    else
      evtwt *= btagsf;
  }

  double ST(htak4.getHT() + zll.at(0).getPt() + goodMet.at(0).getFullPt()); 

  if (goodAK4Jets.at(0).getPt() > 100 ) {
    h1_["cutflow"] -> Fill(7, presel_wt) ; 
    if (goodAK4Jets.at(1).getPt() > 50){
      h1_["cutflow"] -> Fill(8, presel_wt) ; 
      if ( goodBTaggedAK4Jets.size() > 0 ) { 
        h1_["cutflow"] -> Fill(9, evtwt) ;
        if ( ST > STMin_ ) h1_["cutflow"] -> Fill(10, evtwt) ;  
      }
    }
  } //// Completing the cutflow

  vlq::JetCollection  goodWTaggedJets;
  jetWTaggedmaker(evt, goodWTaggedJets);
  cleanjets(goodWTaggedJets, goodMuons); 
  cleanjets(goodWTaggedJets, goodElectrons); 

  vlq::JetCollection goodHTaggedJets; 
  jetHTaggedmaker(evt, goodHTaggedJets);
  cleanjets(goodHTaggedJets, goodMuons); 
  cleanjets(goodHTaggedJets, goodElectrons); 
  
  //// http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_245_v3.pdf
  //// SF of 0.93+/-0.09 required for top tag WP with mistag rate 1% (no subjet b tag): AN2016-245v3
  vlq::JetCollection goodTopTaggedJets;
  jetTopTaggedmaker(evt, goodTopTaggedJets);
  cleanjets(goodTopTaggedJets, goodMuons); 
  cleanjets(goodTopTaggedJets, goodElectrons); 

  double sjbtagsf(1) ;
  double sjbtagsf_bcUp(1) ; 
  double sjbtagsf_bcDown(1) ; 
  double sjbtagsf_lUp(1) ; 
  double sjbtagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodHTaggedJets) {
      csvs.push_back(jet.getCSVSubjet0()) ; 
      csvs.push_back(jet.getCSVSubjet1()) ; 
      pts.push_back(jet.getPtSubjet0()) ; 
      pts.push_back(jet.getPtSubjet1()) ; 
      etas.push_back(jet.getEtaSubjet0()) ; 
      etas.push_back(jet.getEtaSubjet1()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet0()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet1()) ; 
    }

    sjbtagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetHTaggedmaker.idxsjCSVMin_, sjbtagsf, sjbtagsf_bcUp, sjbtagsf_bcDown, sjbtagsf_lUp, sjbtagsf_lDown) ;

  }

  if ( skim_ ) {

    h1_["ht_preSel_bfFit"] -> Fill(htak4.getHT(), presel_wt);
    //fill jet flavor pt and eta for b-tagging efficiencies 
    if (goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){// not sure if we should also cut on ST
      for (vlq::Jet jet : goodAK4Jets) {
        if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
      }
      if ( goodBTaggedAK4Jets.size() > 0 ){
        for (vlq::Jet jet : goodBTaggedAK4Jets) {
          if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
        }
      }
    }

    std::unique_ptr<double> ptr_fevtwt ( new double(evtwt) ) ;
    std::unique_ptr<double> ptr_prewt ( new double(presel_wt) ) ;
    std::unique_ptr<double> ptr_btagsf          ( new double(btagsf         ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcUp     ( new double(btagsf_bcUp    ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcDown   ( new double(btagsf_bcDown  ) ) ;
    std::unique_ptr<double> ptr_btagsf_lUp      ( new double(btagsf_lUp     ) ) ;
    std::unique_ptr<double> ptr_btagsf_lDown    ( new double(btagsf_lDown   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf        ( new double(sjbtagsf       ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcUp   ( new double(sjbtagsf_bcUp  ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcDown ( new double(sjbtagsf_bcDown) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lUp    ( new double(sjbtagsf_lUp   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lDown  ( new double(sjbtagsf_lDown ) ) ;
    evt.put(std::move(ptr_fevtwt), "finalWeight");
    evt.put(std::move(ptr_prewt), "PreWeight");
    evt.put(std::move(ptr_btagsf         ), "btagsf"         );
    evt.put(std::move(ptr_btagsf_bcUp    ), "btagsfbcUp"    );
    evt.put(std::move(ptr_btagsf_bcDown  ), "btagsfbcDown"  );
    evt.put(std::move(ptr_btagsf_lUp     ), "btagsflUp"     );
    evt.put(std::move(ptr_btagsf_lDown   ), "btagsflDown"   );
    evt.put(std::move(ptr_sjbtagsf       ), "sjbtagsf"       );
    evt.put(std::move(ptr_sjbtagsf_bcUp  ), "sjbtagsfbcUp"  );
    evt.put(std::move(ptr_sjbtagsf_bcDown), "sjbtagsfbcDown");
    evt.put(std::move(ptr_sjbtagsf_lUp   ), "sjbtagsflUp"   );
    evt.put(std::move(ptr_sjbtagsf_lDown ), "sjbtagsflDown" );

    if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50 && goodBTaggedAK4Jets.size() > 0 && ST > STMin_){
      std::unique_ptr<vlq::JetCollection> ptr_tjets( new vlq::JetCollection(goodTopTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_wjets( new vlq::JetCollection(goodWTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_hjets( new vlq::JetCollection(goodHTaggedJets) ) ;
      std::unique_ptr<vlq::JetCollection> ptr_bjets( new vlq::JetCollection(goodBTaggedAK4Jets ) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_jets ( new vlq::JetCollection(goodAK4Jets ) ) ; 
      std::unique_ptr<vlq::CandidateCollection> ptr_zllcands ( new vlq::CandidateCollection(zll) ) ; 

      evt.put(std::move(ptr_tjets), "tjets") ; 
      evt.put(std::move(ptr_wjets), "wjets") ; 
      evt.put(std::move(ptr_hjets), "hjets") ; 
      evt.put(std::move(ptr_bjets), "bjets") ; 
      evt.put(std::move(ptr_jets), "jets")  ; 
      evt.put(std::move(ptr_zllcands), "zllcands")  ;
    }   

  } //// if skim 
  else if ( !maketree_ ) { 

    std::string lep("");
    if(zdecayMode_ == "zmumu") {lep = "mu";}
    else if ( zdecayMode_ == "zelel") {lep = "el";}
    else edm::LogError("OS2LAna::filter") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;

    for (auto izll : zll) {
      h1_["mass_z"+lep+lep+"_pre"] -> Fill(izll.getMass(), presel_wt) ;
      h1_["mass_Z"+lep+lep+"_pre"] -> Fill(izll.getMass(), evtwt) ;
      h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ; 
    }
    h1_["nak4_pre"] -> Fill(goodAK4Jets.size(), presel_wt) ;
    h1_["ht_pre"] -> Fill(htak4.getHT(), presel_wt);
    h1_["st_pre"] -> Fill(ST, presel_wt) ;   

    //// Ak4 jet plots
    for(int j=0; j<3; ++j){
      h1_[Form("ptak4jet%d_pre", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), presel_wt) ; 
      h1_[Form("etaak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), presel_wt) ;
      h1_[Form("cvsak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), presel_wt) ;
      h1_[Form("massak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), presel_wt) ;

    }
    h1_["phi_jet1MET_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), presel_wt);

    //// npv
    h1_["npv_noweight_pre"] -> Fill(npv, *h_evtwtGen.product()); 
    h1_["npv_pre"] -> Fill(npv, presel_wt);

    //// met
    h1_["met_pre"] -> Fill(goodMet.at(0).getFullPt(), presel_wt);
    h1_["met1_pre"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
    h1_["metPhi_pre"] -> Fill(goodMet.at(0).getFullPhi(), presel_wt);

    //// Lepton specfic properties
    if ( zdecayMode_ == "zmumu" ){       
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getEta(), presel_wt) ; 
      } 
      for  (unsigned int imu=0; imu<goodMuons.size(); ++imu){
	h1_["Iso04_mu_pre"]->Fill(goodMuons.at(imu).getIso04(), evtwt);
	h1_["Dz_mu_pre"]->Fill(goodMuons.at(imu).getDz(), evtwt);
	h1_["Dxy_mu_pre"]->Fill(goodMuons.at(imu).getDxy(), evtwt);
	h1_["IsGlobalMuon_mu_pre"]->Fill(goodMuons.at(imu).getIsGlobalMuon(), evtwt);
	h1_["IsPFMuon_mu_pre"]->Fill(goodMuons.at(imu).getIsPFMuon(), evtwt);
	h1_["GlbTrkNormChi2_mu_pre"]->Fill(goodMuons.at(imu).getGlbTrkNormChi2(), evtwt);
	h1_["NumberValidMuonHits_mu_pre"]->Fill(goodMuons.at(imu).getNumberValidMuonHits(), evtwt);
	h1_["NumberMatchedStations_mu_pre"]->Fill(goodMuons.at(imu).getNumberMatchedStations(), evtwt);
	h1_["NumberValidPixelHits_mu_pre"]->Fill(goodMuons.at(imu).getNumberValidPixelHits(), evtwt);
	h1_["NumberTrackerLayers_mu_pre"]->Fill(goodMuons.at(imu).getNumberTrackerLayers(), evtwt);


      }
      
      h1_["dr_mumu_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), presel_wt);
    }
    else if (zdecayMode_ == "zelel" ) {
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]   -> Fill(goodElectrons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodElectrons.at(l).getEta(), presel_wt) ; 
      } 

      h1_["dr_elel_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), presel_wt);
      if(additionalPlots_) fillAdditionalPlots(goodElectrons, presel_wt);
    }

    //// Z pt
    for (auto izll : zll) h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ;

    //========================================================
    // Preselection done, proceeding with control selections
    //========================================================

    if ( goodBTaggedAK4Jets.size() == 0) {
      for (auto izll : zll) {
	h1_["nob_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
      }

      h1_["nob_ht"] ->Fill(htak4.getHT(), evtwt);
      h1_["nob_st"] ->Fill(ST, evtwt);

      if (goodMet.at(0).getFullPt()<60){
	h1_["nbjets_met_0btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	h1_["ht_met_0btagcnt"] ->Fill(htak4.getHT(), evtwt);
	h1_["st_met_0btagcnt"] ->Fill(ST, evtwt);
      }
    }
    else if (goodBTaggedAK4Jets.size() > 0){
      for (auto izll : zll) {
	h1_["b_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
	h1_["b_st"] ->Fill(ST, evtwt);
      }
      h1_["1b_ht"] ->Fill(htak4.getHT(), evtwt);
      h1_["1b_st"] ->Fill(ST, evtwt);
      if (goodMet.at(0).getFullPt()<60){
	h1_["nbjets_met_1btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	h1_["ht_met_1btagcnt"] ->Fill(htak4.getHT(), evtwt);
	h1_["st_met_1btagcnt"] ->Fill(ST, evtwt);
      }
    }
    if (goodMet.at(0).getFullPt()<60){
      h1_["nbjets_met_cnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["lowmet_ht"] ->Fill(htak4.getHT(), evtwt);
      h1_["lowmet_st"] ->Fill(ST, evtwt);
    }


    //fill control plots
    if ( goodBTaggedAK4Jets.size() > 0 && ST < STMaxControl_  ) {
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;  
	h1_["mass_Z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;
        h1_["pt_z"+lep+lep+"_cnt"] -> Fill(izll.getPt(), evtwt) ; 
      }
      h1_["nak4_cnt"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht_cnt"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st_cnt"] -> Fill(ST, evtwt) ;   
      h1_["ht1_cnt"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st1_cnt"] -> Fill(ST, evtwt) ;
      h1_["npv_noweight_cnt"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv_cnt"] -> Fill(npv, evtwt);
      h1_["met_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["met1_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi_cnt"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu_cnt"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel_cnt"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }

      //ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d_cnt", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
	h1_[Form("massak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), evtwt) ;

      }
      h1_["phi_jet1MET_cnt"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

    } //// Control region 


    //fill categorization control plots

  if (categorize_){
    if ( goodBTaggedAK4Jets.size() >=1  && ST < 1000) {
      // HCandidate cut(boosted)                                                                             
      for (unsigned i=0; i<goodAK4Jets.size(); i++) { 
	if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
	  TLorentzVector h;
	  h= goodAK4Jets.at(i).getP4();
	  vlq::Candidate h1(h);
	  Hb.push_back(h1);
	}
      }
      //HCandidates (nonboosted)
      HCandsProducer h;
      h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;
      
      //Zcandidate cut (boosted)                                     
      for (unsigned i=0; i<goodAK4Jets.size(); i++) {                                                  
	if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	  TLorentzVector zb;    
	  zb= goodAK4Jets.at(i).getP4();
	  vlq::Candidate zb1(zb);       
	  ZB.push_back(zb1);
	    
	}                                                                                                                                       
      }

      //Z candidate cut (non boosted)                                                                                                            
      ZCandsProducer z;
      z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;

      //Top Candidates (Category D)    
      for (unsigned i=0; i<goodAK4Jets.size(); i++) {
	if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	  TLorentzVector d;
	  d= goodAK4Jets.at(i).getP4();
	  vlq::Candidate d1(d);
	  D.push_back(d1);
	}
      }

      TopCandsProducer top,w;
      // Category BC
      w.operator()(goodAK4Jets, W,B) ;
      for (unsigned i=0; i<W.size(); i++) {
	for (unsigned j=0; j<B.size(); j++) {
	  TLorentzVector bc1;
	  bc1= W.at(i).getP4()+B.at(j).getP4();
	  if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	    vlq::Candidate bc2(bc1);
	    BC.push_back(bc2); 
	  }  
	}
      }
      // cout << "no of good ak4 befor cnt " << goodAK4Jets.size()<<endl;
      // category A
      top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;

      //nb >= 1 and ST <1000 control region
      //  cout << " set 0 " << endl;
      // cout << "no of good ak4 " << goodAK4Jets.size()<<endl;
      // cout << " Hb size is = " << Hb.size()<<endl;
      // cout << "  H size is  = " << H.size()<<endl;
      // cout << " ZB size is = " << ZB.size()<<endl;
      // cout << " Z size is = " << Z.size()<<endl;
      // cout << " D size is = " << D.size()<<endl;
      // cout << " BC size is = " << BC.size()<<endl;
      // cout << " tops size is = " << tops.size()<<endl;
      


      // HCandidate cut(boosted)

      for (unsigned i=0; i<Hb.size(); i++) {
	h1_["H_mass_b_cnt"] -> Fill(Hb.at(i).getMass(), evtwt) ;
	h1_["H_Pt_b_cnt"] -> Fill(Hb.at(i).getPt(), evtwt) ;
      }
      h1_["nHcandidatejets_b_cnt"] -> Fill(Hb.size(), evtwt) ;
      
      // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
      
      //Hcandidate cut(nonboosted)
      for (unsigned i=0; i<H.size(); i++) {
	h1_["H_mass_nb_cnt"] -> Fill(H.at(i).getMass(), evtwt) ;
	h1_["H_Pt_nb_cnt"] -> Fill(H.at(i).getPt(), evtwt) ;
      }
      
      h1_["nHcandidatejets_nb_cnt"] -> Fill(H.size(), evtwt) ;

      //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
      
      double  nHcandidates=0.0;
      double  nHcandidates1=0.0;
      if(Hb.size()>0 ||H.size()>0){
	//h1_["cutflow"] -> Fill(16, evtwt) ;
	nHcandidates = Hb.size()+H.size();
	
	h1_["nHcandidatejets_cnt"] -> Fill(nHcandidates, evtwt) ;
      
      }

      nHcandidates1 = Hb.size()+H.size();
      h1_["nHcandidatejets1_cnt"] -> Fill(nHcandidates1, evtwt) ;

      //Zcandidate cut (boosted)
      for (unsigned i=0; i<ZB.size(); i++) {     
	h1_["Z_mass_a_cnt"] -> Fill(ZB.at(i).getMass(), evtwt) ;
	h1_["Z_Pt_a_cnt"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
      }                                                                                                                                                                                       
      h1_["nzcandidatejets_a_cnt"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
   
      //if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     
      
      //Z candidate cut (non boosted)                                                                                                                                                         
      for (unsigned i=0; i<Z.size(); i++) {
	h1_["Z_mass_b_cnt"] -> Fill(Z.at(i).getMass(), evtwt) ;
	h1_["Z_Pt_b_cnt"] -> Fill(Z.at(i).getPt(), evtwt) ;
      }

      h1_["nzcandidatejets_b_cnt"] -> Fill(Z.size(), evtwt) ;
      
      // if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
 
      double nzcandidates=0.0;
      double nzcandidates1=0.0;
      if ( ZB.size() || Z.size()){
	//h1_["cutflow"] -> Fill(12, evtwt);
	nzcandidates = ZB.size()+ Z.size();
	h1_["nzcandidatejets_tot_cnt"] -> Fill(nzcandidates, evtwt) ;
	
      }
      nzcandidates1 = ZB.size()+ Z.size();
      h1_["nzcandidatejets1_tot_cnt"] -> Fill(nzcandidates1, evtwt) ;

      
      // Category D    
      for (unsigned i=0; i<D.size(); i++) {
	h1_["top_mass_d_cnt"] -> Fill(D.at(i).getMass(), evtwt) ;
	h1_["top_Pt_d_cnt"] -> Fill(D.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_d_cnt"] -> Fill(D.size(), evtwt) ;
      
      // if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}


      // Category BC
      for (unsigned i=0; i<W.size(); i++) {
	h1_["W_mass_bc_cnt"] -> Fill(W.at(i).getMass(), evtwt) ;
      }
      h1_["nWcandidatejets_bc_cnt"] -> Fill(W.size(), evtwt) ;
      
      for (unsigned i=0; i<B.size(); i++) {
	h1_["lightjet_mass_bc_cnt"] -> Fill(B.at(i).getMass(), evtwt) ; 
      }
      h1_["nlightjetcandidatejets_bc_cnt"] -> Fill(B.size(), evtwt) ;
      
      for (unsigned i=0; i<BC.size(); i++) {
	h1_["top_mass_bc_cnt"] -> Fill(BC.at(i).getMass(), evtwt) ;
	h1_["top_Pt_bc_cnt"] -> Fill(BC.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_bc_cnt"] -> Fill(BC.size(), evtwt) ;
      
      // if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}

      for (unsigned i=0; i<tops.size(); i++) {
	h1_["top_mass_a_cnt"] -> Fill(tops.at(i).getMass(), evtwt) ;
	h1_["top_Pt_a_cnt"] -> Fill(tops.at(i).getPt(), evtwt) ;
      }
      h1_["ntopcandidatejets_a_cnt"] -> Fill(tops.size(), evtwt) ;
      
      // if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
      
      double  ntopcandidates=0.0;
      double  ntopcandidates1=0.0;
      if(D.size() || BC.size() ||tops.size()){
	//   h1_["cutflow"] -> Fill(16, evtwt) ;
	ntopcandidates = D.size()+BC.size()+tops.size();
	
	h1_["ntopcandidatejets_cnt"] -> Fill(ntopcandidates, evtwt) ;  
	
      } 
      ntopcandidates1 = D.size()+BC.size()+tops.size();
      h1_["ntopcandidatejets1_cnt"] -> Fill(ntopcandidates1, evtwt) ;
      
      //cout << " pass 2 " <<endl;
      //Z and top corelations and ST tempelates

      //h1_["st_cnt"] -> Fill(ST,evtwt);
      h1_["cutflow1"] -> Fill(1, evtwt) ;
      if (goodBTaggedAK4Jets.size() == 1){
	h1_["cutflow1"] -> Fill(2, evtwt) ;   
	//h1_["st_cnt1b"] -> Fill(ST,evtwt);
      }
      
      if (goodBTaggedAK4Jets.size() >=2){
	//h1_["st_cnt2b"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(3, evtwt) ;
      }

      //n,Z,H,B
      //(1)
      if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	//h1_["st_cntT1Z1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(4, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(8, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(1, evtwt) ;
	    h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(2, evtwt) ;
	    h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(9, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(3, evtwt) ;
	     h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	    
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(4, evtwt) ;
	    h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}  
      }

      //(2)                                                                                                                                                                                        
      if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
	//h1_["st_cntT0Z1"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(5, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(10, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(5, evtwt) ;
	    h1_["st_cntT0Z1H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(6, evtwt) ;
	    h1_["st_cntT0Z1H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(11, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(7, evtwt) ;
	    h1_["st_cntT0Z1H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(8, evtwt) ;
	    h1_["st_cntT0Z1H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }

      //(3)                                                                                                                                                                                      
      if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
	//h1_["st_cntT1Z0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(6, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(12, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(9, evtwt) ;
	    h1_["st_cntT1Z0H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(10, evtwt) ;
	    h1_["st_cntT1Z0H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(13, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(11, evtwt) ;
	     h1_["st_cntT1Z0H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(12, evtwt) ;
	    h1_["st_cntT1Z0H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }

      //(4)                                                                                                                                                                                        
      if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
	//h1_["st_cntT0Z0"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(7, evtwt) ;
	if (nHcandidates >= 1.0){
	  // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(14, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(13, evtwt) ;
	    h1_["st_cntT0Z0H1b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(14, evtwt) ;
	    h1_["st_cntT0Z0H1b2"] -> Fill(ST, evtwt) ;
	  }
	}
	else if (nHcandidates == 0.0){
	  //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(15, evtwt) ;
	  if( goodBTaggedAK4Jets.size() == 1 ){
	    h1_["cutflow2"] -> Fill(15, evtwt) ;
	    h1_["st_cntT0Z0H0b1"] -> Fill(ST, evtwt) ;
	  }
	  else if( goodBTaggedAK4Jets.size() >= 2 ){
	    h1_["cutflow2"] -> Fill(16, evtwt) ;
	    h1_["st_cntT0Z0H0b2"] -> Fill(ST, evtwt) ;
	  }
	}
      }
    }
  }

  Hb.clear();
  H.clear();
  ZB.clear();
  Z.clear();
  D.clear();
  BC.clear();
  tops.clear();









    if ( goodAK4Jets.at(0).getPt() > 100 
        && goodAK4Jets.at(1).getPt() > 50
        && goodBTaggedAK4Jets.size() > 0
        && ST > STMin_ 
       ) { //// Signal region 

      //// fill all the plots in signal region
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;  
	h1_["mass_Z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;
        h1_["pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ; 
      }

      h1_["nak4"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st"] -> Fill(ST, evtwt) ;   
      h1_["npv_noweight"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv"] -> Fill(npv, evtwt);
      h1_["met"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["met1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //// lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }
      //// ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
	h1_[Form("massak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), evtwt) ;  

    }
      h1_["phi_jet1MET"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

      //// fill the b-tagging plots and efficiency maps
      h1_["nbjets"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["ptbjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
      h1_["etabjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;

      //// fill the additional plots
      h1_["nak8"] -> Fill(goodAK8Jets.size(), evtwt) ;
      h1_["nwjet"] -> Fill(goodWTaggedJets.size(), evtwt) ; 
      h1_["nhjet"] -> Fill(goodHTaggedJets.size(), evtwt) ; 
      h1_["ntjet"] -> Fill(goodTopTaggedJets.size(), evtwt) ; 

      if (goodAK8Jets.size() > 0) {
        h1_["ptak8leading"] -> Fill((goodAK8Jets.at(0)).getPt(), evtwt) ; 
        h1_["etaak8leading"] -> Fill((goodAK8Jets.at(0)).getEta(), evtwt) ;
        //h1_["mak8leading"] -> Fill((goodAK8Jets.at(0)).getMass(), evtwt) ; 
        //h1_["trimmedmak8leading"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(), evtwt) ;
        //h1_["prunedmak8leading"] -> Fill((goodAK8Jets.at(0)).getPrunedMass(), evtwt) ;
        h1_["softdropmak8leading"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(), evtwt) ;
      }
      if (goodAK8Jets.size() > 1) {
        h1_["ptak82nd"] -> Fill((goodAK8Jets.at(1)).getPt(), evtwt) ; 
        h1_["etaak82nd"] -> Fill((goodAK8Jets.at(1)).getEta(), evtwt) ;
        //h1_["mak82nd"] -> Fill((goodAK8Jets.at(1)).getMass(), evtwt) ; 
        //h1_["trimmedmak82nd"] -> Fill((goodAK8Jets.at(1)).getTrimmedMass(), evtwt) ;
        //h1_["prunedmak82nd"] -> Fill((goodAK8Jets.at(1)).getPrunedMass(), evtwt) ;
        h1_["softdropmak82nd"] -> Fill((goodAK8Jets.at(1)).getSoftDropMass(), evtwt) ;
      }

 //Do mass reconstruction TPrime
  if (categorize_){

    // HCandidate cut(boosted)                                                                             
    for (unsigned i=0; i<goodAK4Jets.size(); i++) { 
      if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
	TLorentzVector h;
	h= goodAK4Jets.at(i).getP4();
	vlq::Candidate h1(h);
	Hb.push_back(h1);
      }
    }
    //HCandidates (nonboosted)
    HCandsProducer h;
    h.operator()(goodAK4Jets.size(), 2, goodAK4Jets,H)  ;

    //Zcandidate cut (boosted)                                     
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {                                                  
      if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	TLorentzVector zb;    
	zb= goodAK4Jets.at(i).getP4();
	vlq::Candidate zb1(zb);       
	ZB.push_back(zb1);
	
      }                                                                                                                                       
    }

    //Z candidate cut (non boosted)                                                                                                            
    ZCandsProducer z;
    z.operator()(goodAK4Jets.size(), 2, goodAK4Jets,Z) ;
    
    //Top Candidates (Category D)    
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	TLorentzVector d;
	d= goodAK4Jets.at(i).getP4();
	vlq::Candidate d1(d);
	D.push_back(d1);
      }
    }

    TopCandsProducer top,w;
    // Category BC
    w.operator()(goodAK4Jets, W,B) ;
    for (unsigned i=0; i<W.size(); i++) {
      for (unsigned j=0; j<B.size(); j++) {
	TLorentzVector bc1;
	bc1= W.at(i).getP4()+B.at(j).getP4();
	if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	  vlq::Candidate bc2(bc1);
	  BC.push_back(bc2); 
	}  
      }
    }

    // category A
    top.operator()(goodAK4Jets.size(), 3, goodAK4Jets,tops) ;

    for (unsigned i=0; i<Hb.size(); i++) {
      h1_["H_mass_b_sig"] -> Fill(Hb.at(i).getMass(), evtwt) ;
      h1_["H_Pt_b_sig"] -> Fill(Hb.at(i).getPt(), evtwt) ;
    }
    h1_["nHcandidatejets_b_sig"] -> Fill(Hb.size(), evtwt) ;
    
    // if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}

    //Hcandidate cut(nonboosted)
    for (unsigned i=0; i<H.size(); i++) {
      h1_["H_mass_nb_sig"] -> Fill(H.at(i).getMass(), evtwt) ;
      h1_["H_Pt_nb_sig"] -> Fill(H.at(i).getPt(), evtwt) ;
    }
    h1_["nHcandidatejets_nb_sig"] -> Fill(H.size(), evtwt) ;

    //if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    
    double  nHcandidates=0.0;
    double  nHcandidates1=0.0;
    if(Hb.size()>0 ||H.size()>0){
      //h1_["cutflow"] -> Fill(16, evtwt) ;
      nHcandidates = Hb.size()+H.size();
      h1_["nHcandidatejets_sig"] -> Fill(nHcandidates, evtwt) ;
    }

    nHcandidates1 = Hb.size()+H.size();
    h1_["nHcandidatejets1_sig"] -> Fill(nHcandidates1, evtwt) ;
    
    //Zcandidate cut (boosted)
    for (unsigned i=0; i<ZB.size(); i++) {     
      h1_["Z_mass_a_sig"] -> Fill(ZB.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_a_sig"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
    }
    h1_["nzcandidatejets_a_sig"] -> Fill(ZB.size(), evtwt) ;
    
    //if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     
    
    //Z candidate cut (non boosted)
    for (unsigned i=0; i<Z.size(); i++) {
      h1_["Z_mass_b_sig"] -> Fill(Z.at(i).getMass(), evtwt) ;
      h1_["Z_Pt_b_sig"] -> Fill(Z.at(i).getPt(), evtwt) ;
    }
    h1_["nzcandidatejets_b_sig"] -> Fill(Z.size(), evtwt) ;

    // if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
    
    double nzcandidates=0.0;
    double nzcandidates1=0.0;
    if ( ZB.size() || Z.size()){
      // h1_["cutflow"] -> Fill(12, evtwt);
      nzcandidates = ZB.size()+ Z.size();
      h1_["nzcandidatejets_tot_sig"] -> Fill(nzcandidates, evtwt) ;
    }
    nzcandidates1 = ZB.size()+ Z.size();
    h1_["nzcandidatejets1_tot_sig"] -> Fill(nzcandidates1, evtwt) ;

    
    // Category D    
    for (unsigned i=0; i<D.size(); i++) {
      h1_["top_mass_d_sig"] -> Fill(D.at(i).getMass(), evtwt) ;
      h1_["top_Pt_d_sig"] -> Fill(D.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_d_sig"] -> Fill(D.size(), evtwt) ;
    
    // if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}


    // Category BC
    for (unsigned i=0; i<W.size(); i++) {
      h1_["W_mass_bc_sig"] -> Fill(W.at(i).getMass(), evtwt) ;
    }
    h1_["nWcandidatejets_bc_sig"] -> Fill(W.size(), evtwt) ;
    
    for (unsigned i=0; i<B.size(); i++) {
      h1_["lightjet_mass_bc_sig"] -> Fill(B.at(i).getMass(), evtwt) ; 
    }
    h1_["nlightjetcandidatejets_bc_sig"] -> Fill(B.size(), evtwt) ;
    
    for (unsigned i=0; i<BC.size(); i++) {
      h1_["top_mass_bc_sig"] -> Fill(BC.at(i).getMass(), evtwt) ;
      h1_["top_Pt_bc_sig"] -> Fill(BC.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_bc_sig"] -> Fill(BC.size(), evtwt) ;

    // if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}
      
    for (unsigned i=0; i<tops.size(); i++) {
      h1_["top_mass_a_sig"] -> Fill(tops.at(i).getMass(), evtwt) ;
      h1_["top_Pt_a_sig"] -> Fill(tops.at(i).getPt(), evtwt) ;
    }
    h1_["ntopcandidatejets_a_sig"] -> Fill(tops.size(), evtwt) ;
    
    // if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    
    double  ntopcandidates=0.0;
    double  ntopcandidates1=0.0;
    if(D.size() || BC.size() ||tops.size()){
      //   h1_["cutflow"] -> Fill(16, evtwt) ;
      ntopcandidates = D.size()+BC.size()+tops.size();
      h1_["ntopcandidatejets_sig"] -> Fill(ntopcandidates, evtwt) ;  
    } 
    ntopcandidates1 = D.size()+BC.size()+tops.size();
    h1_["ntopcandidatejets1_sig"] -> Fill(ntopcandidates1, evtwt) ;

    // cout << " st befor rebin " << ST <<endl;
    //Z and top corelations and ST tempelates
    if (ST >2500.00) ST = 2498.00;    
    //cout << " st after rebin ****" << ST <<endl;
    h1_["st_sig"] -> Fill(ST,evtwt);
    h1_["cutflow3"] -> Fill(1, evtwt) ;
    if (goodBTaggedAK4Jets.size() == 1){
      h1_["cutflow3"] -> Fill(2, evtwt) ;   
      h1_["st_sig1b"] -> Fill(ST,evtwt);
    }
    if (goodBTaggedAK4Jets.size() >=2){
      h1_["st_sig2b"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(3, evtwt) ;
    }

    //n,Z,H,B
    //(1)
    if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
      h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(4, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(8, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(1, evtwt) ;
	  h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(2, evtwt) ;
	  h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;
	}
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(9, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(3, evtwt) ;
	  h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(4, evtwt) ;
	  h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }  
    }

    //(2)
    if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
      h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(5, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT0Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(10, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(5, evtwt) ;
	  h1_["st_sigT0Z1H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(6, evtwt) ;
	  h1_["st_sigT0Z1H1b2"] -> Fill(ST, evtwt) ;
	}
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT0Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(11, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(7, evtwt) ;
	  h1_["st_sigT0Z1H0b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(8, evtwt) ;
	  h1_["st_sigT0Z1H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(3)
    if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
      h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(6, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT1Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(12, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(9, evtwt) ;
	  h1_["st_sigT1Z0H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(10, evtwt) ;
	  h1_["st_sigT1Z0H1b2"] -> Fill(ST, evtwt) ;
	}
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT1Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(13, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(11, evtwt) ;
	  h1_["st_sigT1Z0H0b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(12, evtwt) ;
	  h1_["st_sigT1Z0H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }

    //(4)
    if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
      h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
      h1_["cutflow3"] -> Fill(7, evtwt) ;
      if (nHcandidates >= 1.0){
	h1_["st_sigT0Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(14, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(13, evtwt) ;
	  h1_["st_sigT0Z0H1b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(14, evtwt) ;
	  h1_["st_sigT0Z0H1b2"] -> Fill(ST, evtwt) ;
	}
      }
      else if (nHcandidates == 0.0){
	h1_["st_sigT0Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(15, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	  h1_["cutflow4"] -> Fill(15, evtwt) ;
	  h1_["st_sigT0Z0H0b1"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	  h1_["cutflow4"] -> Fill(16, evtwt) ;
	  h1_["st_sigT0Z0H0b2"] -> Fill(ST, evtwt) ;
	}
      }
    }
  }












    } //// Signal region 
    else return false;

  } //// if not skim  and not maketree

  if ( maketree_ ) {

    os2ltree_.clearTreeVectors();

    os2ltree_.t_signalType = signalType ; 

    os2ltree_.ta_npv = npv ; 
    os2ltree_.t_evtInfoEventNumber =  evtno ; 
    os2ltree_.t_evtInfoRunNumber = runno ; 
    os2ltree_.t_evtInfoLumiBlock = lumisec ; 

    for (vlq::Electron e : goodElectrons) {
      os2ltree_.t_elPt        .push_back(e.getPt()) ; 
      os2ltree_.t_elPhi       .push_back(e.getPhi());
      os2ltree_.t_elEta       .push_back(e.getEta());
      os2ltree_.t_elE         .push_back(e.getE());
      os2ltree_.t_elCharge    .push_back(e.getCharge());
      os2ltree_.t_elIso03     .push_back(e.getIso03());
      os2ltree_.t_elecIsLoose .push_back(e.getisLoose()); 
      os2ltree_.t_elecIsMedium.push_back(e.getisMedium()); 
      os2ltree_.t_elecIsTight .push_back(e.getisTight()); 
      os2ltree_.t_elecIsVeto  .push_back(e.getisVeto());
    }

    for (vlq::Muon m : goodMuons) {
      os2ltree_.t_muPt         .push_back(m.getPt());
      os2ltree_.t_muPhi        .push_back(m.getPhi());
      os2ltree_.t_muEta        .push_back(m.getEta());
      os2ltree_.t_muE          .push_back(m.getE());
      os2ltree_.t_muCharge     .push_back(m.getCharge());
      os2ltree_.t_muIso04      .push_back(m.getIso04());
      os2ltree_.t_muonIsTight  .push_back(m.getIsTightMuon()); 
      os2ltree_.t_muonIsLoose  .push_back(m.getIsLooseMuon()); 
      os2ltree_.t_muonIsGlobal .push_back(m.getIsGlobalMuon()); 
      os2ltree_.t_muonIsPFMuon .push_back(m.getIsPFMuon()); 
      os2ltree_.t_muonIsTracker.push_back(m.getIsTrackerMuon()); 
    }

    for (vlq::Jet jet : goodAK4Jets) {
      os2ltree_.t_jetAK4Pt           .push_back(jet.getPt());
      os2ltree_.t_jetAK4Phi          .push_back(jet.getPhi());
      os2ltree_.t_jetAK4Eta          .push_back(jet.getEta());
      os2ltree_.t_jetAK4E            .push_back((jet.getP4()).E());
      os2ltree_.t_jetAK4CSV          .push_back(jet.getCSV());
      os2ltree_.t_jetAK4Mass         .push_back(jet.getMass());
      os2ltree_.t_jetAK4HadronFlavour.push_back(jet.getHadronFlavour());
      os2ltree_.t_jetAK4PartonFlavour.push_back(jet.getPartonFlavour());
    }

    for (vlq::Jet bjet : goodBTaggedAK4Jets) {
      os2ltree_.t_jetAK4BPt            .push_back(bjet.getPt());
      os2ltree_.t_jetAK4BPhi           .push_back(bjet.getPhi());
      os2ltree_.t_jetAK4BEta           .push_back(bjet.getEta());
      os2ltree_.t_jetAK4BE             .push_back((bjet.getP4()).E());
      os2ltree_.t_jetAK4BCSV           .push_back(bjet.getCSV());
      os2ltree_.t_jetAK4BMass          .push_back(bjet.getMass());
      os2ltree_.t_jetAK4BHadronFlavour .push_back(bjet.getHadronFlavour());
      os2ltree_.t_jetAK4BPartonFlavour .push_back(bjet.getPartonFlavour());
    }

    for (vlq::Jet ak8 : goodAK8Jets) {
      os2ltree_.t_jetAK8Pt            .push_back(ak8.getPt());
      os2ltree_.t_jetAK8Phi           .push_back(ak8.getPhi());
      os2ltree_.t_jetAK8Eta           .push_back(ak8.getEta());
      os2ltree_.t_jetAK8E             .push_back((ak8.getP4()).E());
      os2ltree_.t_jetAK8CSV           .push_back(ak8.getCSV());
      os2ltree_.t_jetAK8Mass          .push_back(ak8.getMass());
      os2ltree_.t_jetAK8HadronFlavour .push_back(ak8.getHadronFlavour());
      os2ltree_.t_jetAK8PartonFlavour .push_back(ak8.getPartonFlavour());
      os2ltree_.t_jetAK8_tau1         .push_back(ak8.getTau1());
      os2ltree_.t_jetAK8_tau2         .push_back(ak8.getTau2());
      os2ltree_.t_jetAK8_tau3         .push_back(ak8.getTau3());
      os2ltree_.t_jetAK8_MassPruned   .push_back(ak8.getPrunedMass());
      os2ltree_.t_jetAK8_SoftDropMass .push_back(ak8.getSoftDropMass());
      os2ltree_.t_jetAK8_NSubJets     .push_back(ak8.getNSubjets());
    }

    for (vlq::Jet wjet : goodWTaggedJets) {
      os2ltree_.t_jetWJetPt            .push_back(wjet.getPt());
      os2ltree_.t_jetWJetPhi           .push_back(wjet.getPhi());
      os2ltree_.t_jetWJetEta           .push_back(wjet.getEta());
      os2ltree_.t_jetWJetE             .push_back((wjet.getP4()).E());
      os2ltree_.t_jetWJetCSV           .push_back(wjet.getCSV());
      os2ltree_.t_jetWJetMass          .push_back(wjet.getMass());
      os2ltree_.t_jetWJetHadronFlavour .push_back(wjet.getHadronFlavour());
      os2ltree_.t_jetWJetPartonFlavour .push_back(wjet.getPartonFlavour());
      os2ltree_.t_jetWJet_tau1         .push_back(wjet.getTau1());
      os2ltree_.t_jetWJet_tau2         .push_back(wjet.getTau2());
      os2ltree_.t_jetWJet_tau3         .push_back(wjet.getTau3());
      os2ltree_.t_jetWJet_MassPruned   .push_back(wjet.getPrunedMass());
      os2ltree_.t_jetWJet_SoftDropMass .push_back(wjet.getSoftDropMass());
      os2ltree_.t_jetWJet_NSubJets     .push_back(wjet.getNSubjets());
    }

    for (vlq::Jet hjet : goodHTaggedJets) {
      os2ltree_.t_jetHJetPt            .push_back(hjet.getPt());
      os2ltree_.t_jetHJetPhi           .push_back(hjet.getPhi());
      os2ltree_.t_jetHJetEta           .push_back(hjet.getEta());
      os2ltree_.t_jetHJetE             .push_back((hjet.getP4()).E());
      os2ltree_.t_jetHJetCSV           .push_back(hjet.getCSV());
      os2ltree_.t_jetHJetMass          .push_back(hjet.getMass());
      os2ltree_.t_jetHJetHadronFlavour .push_back(hjet.getHadronFlavour());
      os2ltree_.t_jetHJetPartonFlavour .push_back(hjet.getPartonFlavour());
      os2ltree_.t_jetHJet_tau1         .push_back(hjet.getTau1());
      os2ltree_.t_jetHJet_tau2         .push_back(hjet.getTau2());
      os2ltree_.t_jetHJet_tau3         .push_back(hjet.getTau3());
      os2ltree_.t_jetHJet_MassPruned   .push_back(hjet.getPrunedMass());
      os2ltree_.t_jetHJet_SoftDropMass .push_back(hjet.getSoftDropMass());
      os2ltree_.t_jetHJet_NSubJets     .push_back(hjet.getNSubjets());
    }

    for (vlq::Jet tjet : goodTopTaggedJets) {
      os2ltree_.t_jetTopJetPt            .push_back(tjet.getPt());
      os2ltree_.t_jetTopJetPhi           .push_back(tjet.getPhi());
      os2ltree_.t_jetTopJetEta           .push_back(tjet.getEta());
      os2ltree_.t_jetTopJetE             .push_back((tjet.getP4()).E());
      os2ltree_.t_jetTopJetCSV           .push_back(tjet.getCSV());
      os2ltree_.t_jetTopJetMass          .push_back(tjet.getMass());
      os2ltree_.t_jetTopJetHadronFlavour .push_back(tjet.getHadronFlavour());
      os2ltree_.t_jetTopJetPartonFlavour .push_back(tjet.getPartonFlavour());
      os2ltree_.t_jetTopJet_tau1         .push_back(tjet.getTau1());
      os2ltree_.t_jetTopJet_tau2         .push_back(tjet.getTau2());
      os2ltree_.t_jetTopJet_tau3         .push_back(tjet.getTau3());
      os2ltree_.t_jetTopJet_MassPruned   .push_back(tjet.getPrunedMass());
      os2ltree_.t_jetTopJet_SoftDropMass .push_back(tjet.getSoftDropMass());
      os2ltree_.t_jetTopJet_NSubJets     .push_back(tjet.getNSubjets());
    }

    os2ltree_.t_HT = htak4.getHT();
    os2ltree_.t_ST = ST;

    os2ltree_.t_presel_wt = presel_wt;
    os2ltree_.t_evtwt = evtwt;
    os2ltree_.t_btagsf = btagsf;
    os2ltree_.t_btagsf_bcUp = btagsf_bcUp;
    os2ltree_.t_btagsf_bcDown = btagsf_bcDown;
    os2ltree_.t_btagsf_lUp = btagsf_lUp;
    os2ltree_.t_btagsf_lDown = btagsf_lDown;
    os2ltree_.t_sjbtagsf = sjbtagsf;
    os2ltree_.t_sjbtagsf_bcUp = sjbtagsf_bcUp;
    os2ltree_.t_sjbtagsf_bcDown = sjbtagsf_bcDown;
    os2ltree_.t_sjbtagsf_lUp = sjbtagsf_lUp;
    os2ltree_.t_sjbtagsf_lDown = sjbtagsf_lDown;

    if ( evtno < 0) {
      vlq::GenParticleCollection vlqGen = genpart(evt) ;
      for(vlq::GenParticle p : vlqGen) { 
        os2ltree_.t_genPartPt        .push_back(p.getP4().Pt());
        os2ltree_.t_genPartPhi       .push_back(p.getP4().Phi());
        os2ltree_.t_genPartEta       .push_back(p.getP4().Eta());
        os2ltree_.t_genPartE         .push_back(p.getP4().E());
        os2ltree_.t_genPartID        .push_back(p.getPdgID());
        os2ltree_.t_genPartStatus    .push_back(p.getStatus());
        os2ltree_.t_genPartMom1ID    .push_back(p.getMom0PdgID());
        os2ltree_.t_genPartMom2ID    .push_back(p.getMom1PdgID());
        os2ltree_.t_genPartDau1ID    .push_back(p.getDau0PdgID());
        os2ltree_.t_genPartDau2ID    .push_back(p.getDau1PdgID());
        os2ltree_.t_genPartMom1Status.push_back(p.getMom0Status());
        os2ltree_.t_genPartMom2Status.push_back(p.getMom1Status());
        os2ltree_.t_genPartDau1Status.push_back(p.getDau0Status());
        os2ltree_.t_genPartDau2Status.push_back(p.getDau1Status());
      }
    }

    for( vlq::Met m : goodMet ) { 
      os2ltree_.t_metPt .push_back(m.getPt());
      os2ltree_.t_metPhi.push_back(m.getPhi());
      os2ltree_.t_metEta.push_back(m.getEta());
      os2ltree_.t_metE  .push_back(m.getP4().E());
    }

    tree_->Fill(); 
  } //// maketree

  return true ; 
}

// ------------ method called once each job just before starting event loop  ------------
void OS2LAna::beginJob() {
  std::string lep(""); 
  if(zdecayMode_ == "zmumu") {lep = "mu";}
  else if ( zdecayMode_ == "zelel") {lep = "el";}
  else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;


  if (filterSignal_){
    if(skim_ || maketree_){
      const int nCh = 12;
      const char *channel[nCh] = {"bZbZ", "bZbH", "bZtW", "bHbH", "bHtW", "tWtW",
        "tZtZ", "tZtH", "tZbW", "tHtH", "tHbW", "bWbW"};
      h1_["signalEvts_all"] = fs->make<TH1D>("signalEvts_all", "All signal events", 12, 0.5, 12.5) ;
      for (int i=1;i<=nCh;i++) h1_["signalEvts_all"]->GetXaxis()->SetBinLabel(i,channel[i-1]);
    }
    else{
      h1_["signalEvts"] = fs->make<TH1D>("signalEvts", "signal events", 2, 0.5, 2.5) ;
    }
  }

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 10, 0.5, 10.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1 , "All") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2  , "Trig.") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3  , "l^{+}l^{-}") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4 , "75 #lt M(l^{+}l^{-}) #lt 105") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5 , "H_{T} #geq 200") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6 , "N(AK4) #geq 3") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7 , "leading jet pt > 100") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8 , "2nd jet pt > 50") ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9 , "N(b jet) #geq 1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "S_{T} #geq 1000") ; 

  if(skim_){
    // b-tagging efficiency maps:
    h2_["pt_eta_b_all"] = fs->make<TH2D>("pt_eta_b_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_all"] = fs->make<TH2D>("pt_eta_c_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_all"] = fs->make<TH2D>("pt_eta_l_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 

    h2_["pt_eta_b_btagged"] = fs->make<TH2D>("pt_eta_b_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_btagged"] = fs->make<TH2D>("pt_eta_c_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_btagged"] = fs->make<TH2D>("pt_eta_l_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ;

    h1_["ht_preSel_bfFit"] = fs->make<TH1D>("ht_preSel_bfFit", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ; 
  }

  if (!skim_ and !maketree_) {
    TFileDirectory pre = fs->mkdir ("pre");
    TFileDirectory sig = fs->mkdir ("sig");
    TFileDirectory cnt = fs->mkdir ("cnt");
    TFileDirectory cat = fs->mkdir ("cat");
    TFileDirectory *bookDir[4]; bookDir[0] = &pre; bookDir[1] = &cnt; bookDir[2] = &sig;  bookDir[3] = &cat;
    std::vector<string> suffix = {"_pre", "_cnt", ""};

    for (int i=0; i<3; i++){
      h1_[("npv_noweight"+suffix[i]).c_str()] = bookDir[i]->make<TH1D>( ("npv_noweight"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("npv"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("npv"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("nak4"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4"+suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
      h1_[("ht"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("ht"+suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
      h1_[("st"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("st"+suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
      h1_[("met"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met"+suffix[i]).c_str(), ";MET [GeV]", 100, 0., 1000.);
      h1_[("met1"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met1"+suffix[i]).c_str(), ";MET [GeV]", 100, 0., 200.);
      h1_[("metPhi"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("metPhi"+suffix[i]).c_str(), ";#Phi(MET)", 20, -5., 5.);

      //jets
      for(int j=1; j<4; ++j){
        string jetPtName = Form("ptak4jet%d", j)+suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
        h1_[jetPtName.c_str()] = bookDir[i]->make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
        string jetEtaName = Form("etaak4jet%d", j)+suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
        h1_[jetEtaName.c_str()] = bookDir[i]->make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
        string jetCVSName = Form("cvsak4jet%d", j)+suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
        h1_[jetCVSName.c_str()] = bookDir[i]->make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
	string jetMassName = Form("massak4jet%d", j)+suffix[i]; string jetMassTitle  = Form(";Mass(%d leading AK4 jet) ;;",j);
        h1_[jetMassName.c_str()] = bookDir[i]->make<TH1D>(jetMassName.c_str(), jetMassTitle.c_str(), 100 ,0. ,1200.) ;


      }
      string jet1METPhiName = "phi_jet1MET"+suffix[i];
      h1_[jet1METPhiName.c_str()] = bookDir[i]->make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

      //leptons
      // std::string lep("");
      // if(zdecayMode_ == "zmumu") {lep = "mu";}
      // else if ( zdecayMode_ == "zelel") {lep = "el";}
      // else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;
      string mass_Z = "mass_z"+lep+lep+suffix[i];
      h1_[mass_Z.c_str()] = bookDir[i]->make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
      string mass_Z1 = "mass_Z"+lep+lep+suffix[i];
      h1_[mass_Z1.c_str()] = bookDir[i]->make<TH1D>(mass_Z1.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 200, 60., 110.) ;
      string dr_ll = "dr_"+lep+lep+suffix[i];  
      h1_[dr_ll.c_str()] = bookDir[i]->make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      string pt_Z = "pt_z"+lep+lep+suffix[i];
      h1_[pt_Z.c_str()] = bookDir[i]->make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ; 
      for(int l=1; l<3; ++l){
        string lepPtName = "pt_"+lep+Form("%d",l)+suffix[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
        h1_[lepPtName.c_str()] = bookDir[i]->make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
        string lepEtaName = "eta_"+lep+Form("%d",l)+suffix[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
        h1_[lepEtaName.c_str()] = bookDir[i]->make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }
    }

    h1_["ht1_cnt"]   =  cnt.make<TH1D>("ht1_cnt", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    h1_["st1_cnt"]   =  cnt.make<TH1D>("st1_cnt",";S_{T} [GeV]", 200, 0., 1000.) ;

    //additional plots
    h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 
    h1_["ptbjetleading"]  = sig.make<TH1D>("ptbjetleading", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etabjetleading"] = sig.make<TH1D>("etabjetleading", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;

    h1_["nak8"] = sig.make<TH1D>("nak8", ";N(AK8 jets);;" , 11, -0.5,10.5) ; 
    h1_["nwjet"] = sig.make<TH1D>("nwjet", ";N(W jets );;" , 6, -0.5,5.5) ; 
    h1_["nhjet"] = sig.make<TH1D>("nhjet", ";N(H jets );;" , 6, -0.5,5.5) ; 
    h1_["ntjet"] = sig.make<TH1D>("ntjet", ";N(top jets);;" , 6, -0.5,5.5) ; 

    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;


  if(categorize_){

    //for Z and H categories seperetely along with b
    h1_["cutflow1"] = cat.make<TH1D>("cutflow1", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(4, "t1Z1 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(5, "t0Z1") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(6, "t1Z0 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(7, "t0Z0") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(8, "t1Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(9, "t1Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(10, "t0Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(11, "t0Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(13, "t1Z0H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0 " ) ;
      
    h1_["cutflow2"] = cat.make<TH1D>("cutflow2", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;
      
    //for ZH(combined category)
    h1_["cutflow3"] = cat.make<TH1D>("cutflow3", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(4, "t1ZH1 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(5, "t0ZH1") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(6, "t1ZH0 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(7, "t0ZH0") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(8, "t1ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(9, "t1ZH1b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(10, "t0ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(11, "t0ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(12, "t1ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(13, "t1ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(14, "t0ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(15, "t0ZH0b2 " ) ;
      
      
    h1_["cutflow4"] = cat.make<TH1D>("cutflow4", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;
      
    //H candidates                                                                                                           
    h1_["H_mass_b_cnt"]  = cat.make<TH1D>("Hmass-boosted-cnt", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_cnt"]  = cat.make<TH1D>("HPt-boosted-cnt", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_cnt"] = cat.make<TH1D>("nHcandidate-boosted-cnt", ";N(H jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["H_mass_nb_cnt"]  = cat.make<TH1D>("Hmassnb-cnt", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_cnt"]  = cat.make<TH1D>("HPtnb-cnt", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_cnt"] = cat.make<TH1D>("nHcandidatesnb-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
      
    h1_["nHcandidatejets_cnt"] = cat.make<TH1D>("nHcandidates-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_cnt"] = cat.make<TH1D>("nHcandidates1-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
      
    // Z candidates
    h1_["Z_mass_a_cnt"]  = cat.make<TH1D>("Zmass-boosted-cnt", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_cnt"]  = cat.make<TH1D>("ZPt-boosted-cnt", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_cnt"] = cat.make<TH1D>("nzcandidate-boosted-cnt", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["Z_mass_b_cnt"]  = cat.make<TH1D>("Zmass-cnt", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_cnt"]  = cat.make<TH1D>("ZPt-cnt", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_cnt"] = cat.make<TH1D>("nzcandidates-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    h1_["nzcandidatejets_tot_cnt"] = cat.make<TH1D>("nzcandidates-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_cnt"] = cat.make<TH1D>("nzcandidates1-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);  
    // cat A
    h1_["top_mass_a_cnt"]  = cat.make<TH1D>("topmas-A-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_cnt"]  = cat.make<TH1D>("topPt-A-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_cnt"] = cat.make<TH1D>("ntopcandidate-A-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    h1_["top_mass_bc_cnt"]  = cat.make<TH1D>("topmass-Bc-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_cnt"]  = cat.make<TH1D>("topPt-BC-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_cnt"] = cat.make<TH1D>("ntopcandidate-BC-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    // cat D
    h1_["top_mass_d_cnt"]  = cat.make<TH1D>("topmass-D-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_cnt"]  = cat.make<TH1D>("topPt-D-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_cnt"] = cat.make<TH1D>("ntopcandidate-D-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //W and light jet(BC)
    h1_["W_mass_bc_cnt"]  = cat.make<TH1D>("Wmass-BC-cnt", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_cnt"] = cat.make<TH1D>("nWcandidate-BC-cnt", ";N(W candidate jets);;" , 21, -0.5, 20.5);
      
    h1_["lightjet_mass_bc_cnt"]  = cat.make<TH1D>("lightjetmass-BC-cnt", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_cnt"] = cat.make<TH1D>("nlightjetcandidate-cnt", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)
    h1_["ntopcandidatejets_cnt"] = cat.make<TH1D>("ntopcandidate-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_cnt"] = cat.make<TH1D>("ntopcandidate1-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //signal region
      
    //H candidates                                                                                                                                                              
    h1_["H_mass_b_sig"]  = cat.make<TH1D>("Hmass-boosted-sig", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_sig"]  = cat.make<TH1D>("HPt-boosted-sig", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_sig"] = cat.make<TH1D>("nHcandidate-boosted-sig", ";N(H jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["H_mass_nb_sig"]  = cat.make<TH1D>("Hmassnb-sig", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_sig"]  = cat.make<TH1D>("HPtnb-sig", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_sig"] = cat.make<TH1D>("nHcandidatesnb-sig", ";N(H jets);;" , 21, -0.5, 20.5);
      
    h1_["nHcandidatejets_sig"] = cat.make<TH1D>("nHcandidates-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_sig"] = cat.make<TH1D>("nHcandidates1-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);
      
   
    // Z candidates                                                                                                                                                             
    h1_["Z_mass_a_sig"]  = cat.make<TH1D>("Zmass-boosted-sig", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_sig"]  = cat.make<TH1D>("ZPt-boosted-sig", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_sig"] = cat.make<TH1D>("nzcandidate-boosted-sig", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["Z_mass_b_sig"]  = cat.make<TH1D>("Zmass-sig", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_sig"]  = cat.make<TH1D>("ZPt-sig", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_sig"] = cat.make<TH1D>("nzcandidates-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    h1_["nzcandidatejets_tot_sig"] = cat.make<TH1D>("nzcandidates-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_sig"] = cat.make<TH1D>("nzcandidates1-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    // cat A                                                                                                                                                                    
    h1_["top_mass_a_sig"]  = cat.make<TH1D>("topmas-A-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_sig"]  = cat.make<TH1D>("topPt-A-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_sig"] = cat.make<TH1D>("ntopcandidate-A-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    h1_["top_mass_bc_sig"]  = cat.make<TH1D>("topmass-Bc-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_sig"]  = cat.make<TH1D>("topPt-BC-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_sig"] = cat.make<TH1D>("ntopcandidate-BC-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    // cat D                                                                                                                                                                    
    h1_["top_mass_d_sig"]  = cat.make<TH1D>("topmass-D-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_sig"]  = cat.make<TH1D>("topPt-D-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_sig"] = cat.make<TH1D>("ntopcandidate-D-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //W and light jet(BC)                                                                                                                                                       
    h1_["W_mass_bc_sig"]  = cat.make<TH1D>("Wmass-BC-sig", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_sig"] = cat.make<TH1D>("nWcandidate-BC-sig", ";N(W candidate jets);;" , 21, -0.5, 20.5);
      
    h1_["lightjet_mass_bc_sig"]  = cat.make<TH1D>("lightjetmass-BC-sig", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_sig"] = cat.make<TH1D>("nlightjetcandidate-sig", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)                                                                                                                                                       
    h1_["ntopcandidatejets_sig"] = cat.make<TH1D>("ntopcandidate-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_sig"] = cat.make<TH1D>("ntopcandidate1-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);
   
    //ST plots
      
    //top,Z, Higgs
      
    h1_["st_sig"] =cat.make<TH1D>("ST_sig", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sig1b"] =cat.make<TH1D>("ST_sig1b", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sig2b"] =cat.make<TH1D>("ST_sig2b", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT1Z1"] =cat.make<TH1D>("ST_sigT1Z1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1"] =cat.make<TH1D>("ST_sigT0Z1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0"] =cat.make<TH1D>("ST_sigT1Z0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0"] =cat.make<TH1D>("ST_sigT0Z0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT1Z1H1"] =cat.make<TH1D>("ST_sigT1Z1H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0"] =cat.make<TH1D>("ST_sigT1Z1H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1"] =cat.make<TH1D>("ST_sigT0Z1H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0"] =cat.make<TH1D>("ST_sigT0Z1H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1"] =cat.make<TH1D>("ST_sigT1Z0H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0"] =cat.make<TH1D>("ST_sigT1Z0H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1"] =cat.make<TH1D>("ST_sigT0Z0H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0"] =cat.make<TH1D>("ST_sigT0Z0H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
      
    h1_["st_sigT1Z1H1b1"] =cat.make<TH1D>("ST_sigT1Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H1b2"] =cat.make<TH1D>("ST_sigT1Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b1"] =cat.make<TH1D>("ST_sigT1Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b2"] =cat.make<TH1D>("ST_sigT1Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT0Z1H1b1"] =cat.make<TH1D>("ST_sigT0Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1b2"] =cat.make<TH1D>("ST_sigT0Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b1"] =cat.make<TH1D>("ST_sigT0Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b2"] =cat.make<TH1D>("ST_sigT0Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT1Z0H1b1"] =cat.make<TH1D>("ST_sigT1Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1b2"] =cat.make<TH1D>("ST_sigT1Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b1"] =cat.make<TH1D>("ST_sigT1Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b2"] =cat.make<TH1D>("ST_sigT1Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT0Z0H1b1"] =cat.make<TH1D>("ST_sigT0Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1b2"] =cat.make<TH1D>("ST_sigT0Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b1"] =cat.make<TH1D>("ST_sigT0Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b2"] =cat.make<TH1D>("ST_sigT0Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_cntT1Z1H1b1"] =cat.make<TH1D>("ST_cntT1Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z1H1b2"] =cat.make<TH1D>("ST_cntT1Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z1H0b1"] =cat.make<TH1D>("ST_cntT1Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z1H0b2"] =cat.make<TH1D>("ST_cntT1Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_cntT0Z1H1b1"] =cat.make<TH1D>("ST_cntT0Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z1H1b2"] =cat.make<TH1D>("ST_cntT0Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z1H0b1"] =cat.make<TH1D>("ST_cntT0Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z1H0b2"] =cat.make<TH1D>("ST_cntT0Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_cntT1Z0H1b1"] =cat.make<TH1D>("ST_cntT1Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z0H1b2"] =cat.make<TH1D>("ST_cntT1Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z0H0b1"] =cat.make<TH1D>("ST_cntT1Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT1Z0H0b2"] =cat.make<TH1D>("ST_cntT1Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_cntT0Z0H1b1"] =cat.make<TH1D>("ST_cntT0Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z0H1b2"] =cat.make<TH1D>("ST_cntT0Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z0H0b1"] =cat.make<TH1D>("ST_cntT0Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_cntT0Z0H0b2"] =cat.make<TH1D>("ST_cntT0Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.); 





 }

  //additional plots
  //o batg region plots
  h1_["nob_ht"]= cnt.make<TH1D>("nob_ht", "HT no b tags", 100, 0., 3000.);
  h1_["nob_st"] = fs->make<TH1D>("nob_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["nob_pt_z"+lep+lep] = fs->make<TH1D>(("nob_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_pt_z"+lep+lep] = cnt.make<TH1D>(("b_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_st"] = cnt.make<TH1D>("b_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["pt_zlight_pre"] = fs->make<TH1D>("pt_zlight_pre", "p_{T} (Z + q_{light}) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zb_pre"] = fs->make<TH1D>("pt_zb_pre", "p_{T} (Z + b) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zc_pre"] = fs->make<TH1D>("pt_zc_pre", "p_{T} (Z + c) [GeV]", 100, 0., 2000.) ;
    





    ////electrons specific varaibles in EE and EB at preselection level
    if (zdecayMode_ == "zelel" && additionalPlots_){
      h1_["Eta_EB_el_pre"] = pre.make<TH1D>("Eta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
      h1_["Eta_EE_el_pre"] = pre.make<TH1D>("Eta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ;
      h1_["Iso03_EB_el_pre"] = pre.make<TH1D>("Iso03_EB_el_pre", ";Iso03 (EB);;", 100,0,0.3) ;
      h1_["Iso03_EE_el_pre"] = pre.make<TH1D>("Iso03_EE_el_pre", ";Iso03 (EE);;", 100,0,0.3) ;
      h1_["dEtaInSeed_EB_el_pre"] = pre.make<TH1D>("dEtaInSeed_EB_el_pre", ";dEtaInSeed (EB);;", 200,-0.05,0.05) ;
      h1_["dEtaInSeed_EE_el_pre"] = pre.make<TH1D>("dEtaInSeed_EE_el_pre", ";dEtaInSeed (EE);;", 200,-0.05,0.05) ;
      h1_["dPhiIn_EB_el_pre"] = pre.make<TH1D>("dPhiIn_EB_el_pre", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
      h1_["dPhiIn_EE_el_pre"] = pre.make<TH1D>("dPhiIn_EE_el_pre", ";dPhiIn (EE);;", 100,-0.2,0.2);
      h1_["Dz_EB_el_pre"] = pre.make<TH1D>("Dz_EB_el_pre",";dZ (EB);;", 200,-0.1,0.1) ;
      h1_["Dz_EE_el_pre"] = pre.make<TH1D>("Dz_EE_el_pre", ";dZ (EE);;", 200,-0.4,0.4) ;
      h1_["Dxy_EB_el_pre"] = pre.make<TH1D>("Dxy_EB_el_pre", ";d0 (EB);;", 100,-0.1,0.1) ;
      h1_["Dxy_EE_el_pre"] = pre.make<TH1D>("Dxy_EE_el_pre", ";d0 (EE);;", 100,-0.1,0.1) ;
      h1_["SCETA_EB_el_pre"] = pre.make<TH1D>("scEta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
      h1_["SCETA_EE_el_pre"] = pre.make<TH1D>("scEta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ; 
      h1_["Full5x5siee_EB_el_pre"] = pre.make<TH1D>("Full5x5siee_EB_el_pre", ";Full5X5SigmaIEtaIEta (EB);;", 200,0,0.01) ;
      h1_["Full5x5siee_EE_el_pre"] = pre.make<TH1D>("Full5x5siee_EE_el_pre", ";Full5X5SigmaIEtaIEta (EE);;", 100,0,0.03) ;
      h1_["HoE_EB_el_pre"] = pre.make<TH1D>("HoE_EB_el_pre", ";H/E (EB);;", 200,0,0.05) ;
      h1_["HoE_EE_el_pre"] = pre.make<TH1D>("HoE_EE_el_pre", ";H/E (EE);;", 200,0,0.1) ;
      h1_["ooEmooP_EB_el_pre"] = pre.make<TH1D>("ooEmooP_EB_el_pre", ";(1/E - 1/P) (EB);;", 200,0,0.02) ;
      h1_["ooEmooP_EE_el_pre"] = pre.make<TH1D>("ooEmooP_EE_el_pre", ";(1/E - 1/P) (EE);;", 200,0,0.02) ;
      h1_["missHits_EB_el_pre"] = pre.make<TH1D>("missHits_EB_el_pre", ";Expected missing Hits (EB);;", 4,-0.5,3.5) ;
      h1_["missHits_EE_el_pre"] = pre.make<TH1D>("missHits_EE_el_pre", ";Expected missing Hits (EE);;", 4,-0.5,3.5) ;
      h1_["conveto_EB_el_pre"] = pre.make<TH1D>("conveto_EB_El_pre", ";has matched Conveto (EB);;", 4,-0.5,3.5) ;
      h1_["conveto_EE_el_pre"] = pre.make<TH1D>("conveto_EE_el_pre", ";has matched Conveto (EE);;", 4,-0.5,3.5) ;
    }
 
  if (zdecayMode_ == "zmumu"){

    h1_["Iso04_mu_pre"] = pre.make<TH1D>("Iso04_mu_pre", ";Iso04 ;;", 100,0,0.3) ;
    h1_["Dz_mu_pre"] = pre.make<TH1D>("Dz_mu_pre",";dZ ;;", 200,-0.1,0.1) ;
    h1_["Dxy_mu_pre"] = pre.make<TH1D>("Dxy_mu_pre", ";dxy;;", 100,-0.1,0.1) ;
    h1_["IsGlobalMuon_mu_pre"] = pre.make<TH1D>("IsGlobalMuon_mu_pre", ";IsGlobalMuon;;", 4,-0.5,3.5) ;
    h1_["IsPFMuon_mu_pre"] = pre.make<TH1D>("IsPFMuon_mu_pre", ";IsPFMuon;;", 4,-0.5,3.5) ;
    h1_["GlbTrkNormChi2_mu_pre"] = pre.make<TH1D>("GlbTrkNormChi2_mu_pre", ";GlbTrkNormChi2;;", 200,0.0,10.0) ;
    h1_["NumberValidMuonHits_mu_pre"] = pre.make<TH1D>("NumberValidMuonHits_mu_pre", ";NumberValidMuonHits;;", 20,-0.5,19.5) ;
    h1_["NumberMatchedStations_mu_pre"] = pre.make<TH1D>("NumberMatchedStations_mu_pre", ";NumberMatchedStations;;", 20,-0.5,19.5) ;
    h1_["NumberValidPixelHits_mu_pre"] = pre.make<TH1D>("NumberValidPixelHits_mu_pre", ";NumberValidPixelHits;;", 20,-0.5,19.5) ;
    h1_["NumberTrackerLayers_mu_pre"] = pre.make<TH1D>("NumberTrackerLayers_mu_pre", ";NumberTrackerLayers;;", 20,-0.5,19.5) ;


  }
    //additional plots
    h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 

    //addition nb plots
    h1_["nbjets_met_sig"] = sig.make<TH1D>("nbjets_sig", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_cnt"] = sig.make<TH1D>("nbjets_cnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_0btagcnt"] = sig.make<TH1D>("nbjets_0btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_1btagcnt"] = sig.make<TH1D>("nbjets_1btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
  
    h1_["ht_met_0btagcnt"]   =  sig.make<TH1D>( "ht_0btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["1b_ht"]   =  sig.make<TH1D>( "1b_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["ht_met_1btagcnt"]   =  sig.make<TH1D>( "ht_1btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_ht"]   =  sig.make<TH1D>("lowmet_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    


    h1_["st_met_0btagcnt"]   =  sig.make<TH1D>( "st_0btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["1b_st"]   =  sig.make<TH1D>( "1b_st", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["st_met_1btagcnt"]   =  sig.make<TH1D>( "st_1btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_st"]   =  sig.make<TH1D>("lowmet_st", ";S_{T} [GeV]", 100, 0., 4000.) ;




  }

  if (maketree_) {
    tree_ = fs->make<TTree>("tree", "HH4b") ; 
    os2ltree_.RegisterTree(tree_) ; 
  } 

}

void OS2LAna::endJob() {

  return ; 
}

DEFINE_FWK_MODULE(OS2LAna);
