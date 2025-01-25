//syuusei tochu
// my Analysis include
#include <iostream>
#include "TVector3.h"
#include "TList.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
// #include "PWGDQ/Core/VarManager.h"
// #include "PWGDQ/Core/HistogramManager.h"
// #include "PWGDQ/Core/AnalysisCut.h"
// #include "PWGDQ/Core/AnalysisCompositeCut.h"
// #include "PWGDQ/Core/HistogramsLibrary.h"
// #include "PWGDQ/Core/CutsLibrary.h"
// #include "PWGDQ/Core/MCSignal.h"
// #include "PWGDQ/Core/MCSignalLibrary.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/TableHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMFTTracks = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>;

struct DecayChainTracer {
    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
    // histogram and configurable
    // configurable
    Configurable<int> targetPdgCode{"targetPdgCode", 333, "Target particle PDG code for decay chain tracing"}; // phi meson
    Configurable<int> fHist_pt_width{"fHist_pt_width", 3000, "histogram pt width"};
    Configurable<double> fHist_pt_min{"fHist_pt_min", 0, "histogram pt min"};
    Configurable<double> fHist_pt_max{"fHist_pt_max", 30, "histogram pt max"};

    Configurable<int> fHist_eta_width{"fHist_Eta_width", 800, "histogram eta width"};
    Configurable<double> fHist_eta_min{"fHist_Eta_min", -6, "histogram eta min"};
    Configurable<double> fHist_eta_max{"fHist_Eta_max", -2, "histogram eta max"};

    Configurable<int> fHist_phi_width{"fHist_Phi_width", 80, "histogram phi width"};
    Configurable<double> fHist_phi_min{"fHist_Phi_min", -4, "histogram phi min"};
    Configurable<double> fHist_phi_max{"fHist_Phi_max", 4, "histogram phi max"};

        
    void init(InitContext const&) {
        const AxisSpec axisPDGcode{10000,-5000.5,4999.5};
        const AxisSpec axisPt{fHist_pt_width, fHist_pt_min, fHist_pt_max};
        const AxisSpec axisEta{fHist_eta_width, fHist_eta_min, fHist_eta_max};
        const AxisSpec axisPhi{fHist_phi_width, fHist_phi_min, fHist_phi_max};
        const AxisSpec axisPosZ{40,-20,20};
        const AxisSpec axisNum{10,-0.5,9.5};
        // define histogram for test
        //histos.add("mothers_vs_daughters","mothers_vs_daughters",kTH2D, {axisNum,axisNum});
        //histos.add("mothers","mothers",kTH1D, {axisNum});
        //histos.add("daughters","daughters",kTH1D, {axisNum});
        //histos.add("PDGcode_mother","PDGcode_mother",kTH1D,{axisPDGcode});
        //histos.add("PDGcode_daughters","PDGcode_daughters",kTH1D,{axisPDGcode});
        //histos.add("nomother_PDGcode_daughters","nomother_PDGcode_daughters",kTH1D,{axisPDGcode});
        histos.add("PDGcode_all_mcTrack","PDGcode_all_mcTrack",kTH1D,{axisPDGcode});
        //histos.add("PDGcode_phi_mcTrack","PDGcode_phi_mcTrack",kTH1D,{axisPDGcode});
        // define histogram mcEvent
        histos.add("mcEvent_posZ", "mcEvent_posZ", kTH1D, {axisPosZ});        
        // define histogram mcTrack
        histos.add("mcTrack_Phi_pt", "mcTrack_Phi_pt", kTH1D, {axisPt});
        histos.add("mcTrack_Phi_phi", "mcTrack_Phi_phi", kTH1D, {axisPhi});
        histos.add("mcTrack_Phi_eta", "mcTrack_Phi_eta", kTH1D, {axisEta});

        // define histogram muon
        histos.add("muon_pt", "muon_pt", kTH1D, {axisPt});
        histos.add("muon_phi", "muon_phi", kTH1D, {axisPhi});
        histos.add("muon_eta", "muon_eta", kTH1D, {axisEta});

        histos.add("muon_Phi_pt", "muon_Phi_pt", kTH1D, {axisPt});
        histos.add("muon_Phi_phi", "muon_Phi_phi", kTH1D, {axisPhi});
        histos.add("muon_Phi_eta", "muon_Phi_eta", kTH1D, {axisEta});

        // define histogram for collision
        // define histogram for fwdtrack
    }
    
    void calGenPhi(aod::McCollisions const& mcEvents,
               aod::McParticles const& mcTracks)
    {
        // calculate MC gen track
    
        for (auto mcEvent: mcEvents){
            histos.fill(HIST("mcEvent_posZ"),mcEvent.posZ());
        }

        for (auto mcTrack : mcTracks){
            histos.fill(HIST("PDGcode_all_mcTrack"),mcTrack.pdgCode());
            if(mcTrack.pdgCode() == targetPdgCode){
                if (mcTrack.has_mothers()){
                    // MB 
                    if (mcTrack.has_daughters()) {
                        //histos.fill(HIST("PDGcode_phi_mcTrack"),mcTrack.pdgCode());
                        auto mothers = mcTrack.mothers_as<aod::McParticles>();
                        auto daughters = mcTrack.daughters_as<aod::McParticles>();

                        // check mother and daughter PDG
                        int countP = 0;
                        int countmu = 0;
                        //histos.fill(HIST("mothers"),mothers.size());
                            for(auto mother : mothers){
                                //histos.fill(HIST("PDGcode_mother"),mother.pdgCode());
                                if(mother.pdgCode()==2212){
                                    countP++;
                                } 
                            };
                        
                        //histos.fill(HIST("daughters"),daughters.size());
                            for(auto daughter : daughters){
                                //histos.fill(HIST("PDGcode_daughters"),daughter.pdgCode());
                                //std::cout << "daughterPDG = " << daughter.pdgCode()<< std::endl; 
                                if(abs(daughter.pdgCode())==13){
                                    countmu++;
                                    //std::cout << "daughterPDG = 13 " << std::endl; 
                                }
                            };
                            //std::cout << "==================" << std::endl; 
                        if(countmu == 2){
                            histos.fill(HIST("mcTrack_Phi_pt"), mcTrack.pt());
                            histos.fill(HIST("mcTrack_Phi_phi"), mcTrack.phi());
                            histos.fill(HIST("mcTrack_Phi_eta"), mcTrack.eta());
                        }
                    }// if has_daughter

                }else{
                    // entered particles
                    if(mcTrack.has_daughters()){
                        auto daughters = mcTrack.daughters_as<aod::McParticles>();
                        int countmu_nomother = 0;
                        for(auto daughter : daughters){
                            //histos.fill(HIST("nomother_PDGcode_daughters"),daughter.pdgCode());
                            //std::cout << "daughterPDG = " << daughter.pdgCode()<< std::endl; 
                            if(abs(daughter.pdgCode())==13){
                                countmu_nomother++;
                                //std::cout << "daughterPDG = 13 " << std::endl; 
                            }
                        };
                        //std::cout << "==================" << std::endl; 

                        if(countmu_nomother == 2){
                            histos.fill(HIST("mcTrack_Phi_pt"), mcTrack.pt());
                            histos.fill(HIST("mcTrack_Phi_phi"), mcTrack.phi());
                            histos.fill(HIST("mcTrack_Phi_eta"), mcTrack.eta());
                        }
                    }// end if has daughters

                }// end if has mothers
            }// end if pdgCode == 13
        }// end mcTracks for
        std::cout << "finish calGenPhi" << std::endl;
    }; // end void traceDecay

    template <typename TEvents, typename TMuons>
    void calRecoPhi(TEvents const& events,
                    aod::FwdTrackAssoc const& muonAssocs,
                    TMuons const& muons){

        for(auto assoc : muonAssocs){
            auto muon = assoc.template fwdtrack_as<TMuons>();
            if(muon.has_mcParticle()){
                auto mcMuon = muon.template mcParticle_as<aod::McParticles>();
                auto mothers = mcMuon.template mothers_as<aod::McParticles>();
                int countTarget = 0;
                for(auto mother : mothers){
                    std::cout << "motherPDG = " << mother.pdgCode()<< std::endl;
                    if(mother.pdgCode() == targetPdgCode){
                        countTarget++;
                    }
                }
                std::cout << "+++++++++++++++++++ " << std::endl; 
                // mother <- phi only
                    if(countTarget == 1){
                        histos.fill(HIST("muon_pt"), muon.pt());
                        histos.fill(HIST("muon_phi"), muon.phi());
                        histos.fill(HIST("muon_eta"), muon.eta());
                    }
            };// end has_mcParticles
        }

        std::cout << "finish calRecoPhi" << std::endl;
    }

    // main process
    void process(MyEvents const& events,
               aod::FwdTrackAssoc const& muonAssocs,
               MyMuonsWithCov const& muons,
               aod::McCollisions const& mcEvents,
               aod::McParticles const& mcTracks){
        calGenPhi(mcEvents,mcTracks);
        calRecoPhi(events,muonAssocs,muons);
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& config) {
  return WorkflowSpec{
      adaptAnalysisTask<DecayChainTracer>(config),
  };
}
