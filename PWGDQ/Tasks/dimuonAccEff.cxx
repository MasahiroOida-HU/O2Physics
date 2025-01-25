// my Analysis include
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
// other includes
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "MathUtils/Primitive2D.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "EventFiltering/Zorro.h"

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

Preslice<aod::FwdTrackAssoc> perCollision = aod::track_association::collisionId;




struct DimuonAccEff {
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

    Configurable<int> fHist_mass_width{"fHist_mass_width", 20, "histogram mass width"};
    Configurable<double> fHist_mass_min{"fHist_mass_min", 0, "histogram mass min"};
    Configurable<double> fHist_mass_max{"fHist_mass_max", 2, "histogram mass max"};
        
    void init(InitContext const&) {
        const AxisSpec axisPDGcode{10000,-5000.5,4999.5};
        const AxisSpec axisPt{fHist_pt_width, fHist_pt_min, fHist_pt_max};
        const AxisSpec axisEta{fHist_eta_width, fHist_eta_min, fHist_eta_max};
        const AxisSpec axisPhi{fHist_phi_width, fHist_phi_min, fHist_phi_max};
        const AxisSpec axisPosZ{40,-20,20};
        const AxisSpec axisNum{10,-0.5,9.5};
        const AxisSpec axisMass{fHist_mass_width, fHist_mass_min, fHist_mass_max};
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
        // define histogram for dimuon
        histos.add("dimuon_mass_SEPM", "dimuon_mass_SEPM", kTH1D, {axisMass});
        histos.add("dimuon_mass_SEPP", "dimuon_mass_SEPP", kTH1D, {axisMass});
        histos.add("dimuon_mass_SEMM", "dimuon_mass_SEMM", kTH1D, {axisMass});
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
                    TMuons const& /*muons*/
                    ){

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

        double MuonMass = 0.1056584;//GeV/c
        double dimuonMass = 0;
        //same event pairing
        for (auto event : events){
             auto groupedAssocs = muonAssocs.sliceBy(perCollision, event.globalIndex());
            int sign1 = 0;
            int sign2 = 0;
            for (auto& [a1, a2] : o2::soa::combinations(groupedAssocs, groupedAssocs)) {
                auto t1 = a1.template fwdtrack_as<TMuons>();
                auto t2 = a2.template fwdtrack_as<TMuons>();
                if (t1.matchMCHTrackId() == t2.matchMCHTrackId())
                    continue;
                if (t1.matchMFTTrackId() == t2.matchMFTTrackId())
                    continue;
                sign1 = t1.sign();
                sign2 = t2.sign();

                //calculate dimuon invariant mass
                ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), MuonMass);
                ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), MuonMass);
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                dimuonMass = v12.M();
                
                if(sign1 * sign2 > 0){
                    if(sign1 > 0){
                        histos.fill(HIST("dimuon_mass_SEPP"),dimuonMass);
                    }else if(sign1 < 0){
                        histos.fill(HIST("dimuon_mass_SEMM"),dimuonMass);
                    }
                }else if(sign1 * sign2 < 0){
                    histos.fill(HIST("dimuon_mass_SEPM"),dimuonMass);
                }

                //VarManager::FillPairVertexing(event, t1, t2, true);
            }//end combination
        }//end event

        std::cout << "finish calRecoPhi" << std::endl;
    };

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
      adaptAnalysisTask<DimuonAccEff>(config),
  };
}
