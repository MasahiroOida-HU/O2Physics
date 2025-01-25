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
#include "CCDB/CcdbApi.h"
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
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "MathUtils/Primitive2D.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "EventFiltering/Zorro.h"

#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

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

using BCsWithTimestamps = soa::Join<o2::aod::BCs, o2::aod::Timestamps>;

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

    Configurable<double> fmatchingchi2cut{"fmatchingchi2cut", 3000, "single muon cut: matching chi2"};
    Configurable<double> fMuonMatchEtaMax{"fMuonMatchEtaMax", -2.5, "Refit eta max cut"};
    Configurable<double> fMuonMatchEtaMin{"fMuonMatchEtaMib", -3.6, "Refit eta min cut"};
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
        Service<o2::ccdb::BasicCCDBManager> ccdb;
        o2::ccdb::CcdbApi ccdbApi;
        o2::field::MagneticField* fieldB;
        Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
        Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    int mRunNumber;
        
    void init(InitContext const&) {
        ccdb->setURL(fConfigCcdbUrl);
        ccdb->setCaching(true);
        ccdb->setLocalObjectValidityChecking();
        ccdb->setFatalWhenNull(false);
        ccdbApi.init(fConfigCcdbUrl);
        const AxisSpec axisPDGcode{10000,-5000.5,4999.5};
        const AxisSpec axisPt{fHist_pt_width, fHist_pt_min, fHist_pt_max};
        const AxisSpec axisEta{fHist_eta_width, fHist_eta_min, fHist_eta_max};
        const AxisSpec axisPhi{fHist_phi_width, fHist_phi_min, fHist_phi_max};
        const AxisSpec axisPosZ{40,-20,20};
        const AxisSpec axisNum{10,-0.5,9.5};
        const AxisSpec axisMass{fHist_mass_width, fHist_mass_min, fHist_mass_max};
        const AxisSpec axispDCA{10000,0,1000};
        const AxisSpec axisRabs{1000,0,100};
        const AxisSpec axischi2{300,0,300};
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
        histos.add("mcTrack_Phi_PtvsEta", "mcTrack_Phi_PtvsEta", kTH2D, {axisPt,axisEta});

        // define histogram muon
        histos.add("muon_pt", "muon_pt", kTH1D, {axisPt});
        histos.add("muon_phi", "muon_phi", kTH1D, {axisPhi});
        histos.add("muon_eta", "muon_eta", kTH1D, {axisEta});
        histos.add("muon_mother", "muon_mother", kTH1D, {axisPDGcode});
        
        histos.add("muon_BeforCut_trackType", "muon_BeforCut_trackType", kTH1D, {axisNum});
        histos.add("muon_BeforCut_pt", "muon_BeforCut_pt", kTH1D, {axisPt});
        histos.add("muon_BeforCut_phi", "muon_BeforCut_phi", kTH1D, {axisPhi});
        histos.add("muon_BeforCut_eta", "muon_BeforCut_eta", kTH1D, {axisEta});
        histos.add("muon_BeforCut_pDCA", "muon_BeforCut_pDCA", kTH1D, {axispDCA});
        histos.add("muon_BeforCut_Rabs", "muon_BeforCut_Rabs", kTH1D, {axisRabs});
        histos.add("muon_BeforCut_chi2", "muon_BeforCut_chi2", kTH1D, {axischi2});

        histos.add("muon_AfterCut_trackType", "muon_AfterCut_trackType", kTH1D, {axisNum});
        histos.add("muon_AfterCut_pt", "muon_AfterCut_pt", kTH1D, {axisPt});
        histos.add("muon_AfterCut_phi", "muon_AfterCut_phi", kTH1D, {axisPhi});
        histos.add("muon_AfterCut_eta", "muon_AfterCut_eta", kTH1D, {axisEta});
        histos.add("muon_AfterCut_pDCA", "muon_AfterCut_pDCA", kTH1D, {axispDCA});
        histos.add("muon_AfterCut_Rabs", "muon_AfterCut_Rabs", kTH1D, {axisRabs});
        histos.add("muon_AfterCut_chi2", "muon_AfterCut_chi2", kTH1D, {axischi2});

            //histos.add("muon_Phi_pt", "muon_Phi_pt", kTH1D, {axisPt});
            //histos.add("muon_Phi_phi", "muon_Phi_phi", kTH1D, {axisPhi});
            //histos.add("muon_Phi_eta", "muon_Phi_eta", kTH1D, {axisEta});
        // define histogram for collision
        // define histogram for fwdtrack
        // define histogram for dimuon
        histos.add("dimuon_Phi_PtvsEta","dimuon_Phi_PtvsEta",kTH2D, {axisPt,axisEta});
        histos.add("dimuon_mass_SEPM", "dimuon_mass_SEPM", kTH1D, {axisMass});
        histos.add("dimuon_mass_SEPP", "dimuon_mass_SEPP", kTH1D, {axisMass});
        histos.add("dimuon_mass_SEMM", "dimuon_mass_SEMM", kTH1D, {axisMass});
    }

    void initCCDB(BCsWithTimestamps::iterator const& bc){
        if (mRunNumber == bc.runNumber()) {
        return;
        }
        mRunNumber = bc.runNumber();
        std::map<string, string> metadata;
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
        auto ts = soreor.first;
        auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
        o2::base::Propagator::initFieldFromGRP(grpmag);
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
        }
        o2::mch::TrackExtrap::setField();
        fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
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
                        auto mothers = mcTrack.mothers_as<aod::McParticles>();
                        auto daughters = mcTrack.daughters_as<aod::McParticles>();

                        // check mother and daughter PDG
                        int countP = 0;
                        int countmu = 0;
                            for(auto mother : mothers){
                                if(mother.pdgCode()==2212){
                                    countP++;
                                } 
                            };
                        
                            for(auto daughter : daughters){

                                if(abs(daughter.pdgCode())==13){
                                    countmu++;
                                }
                            };
                        if(countmu == 2){
                            histos.fill(HIST("mcTrack_Phi_pt"), mcTrack.pt());
                            histos.fill(HIST("mcTrack_Phi_phi"), mcTrack.phi());
                            histos.fill(HIST("mcTrack_Phi_eta"), mcTrack.eta());
                            histos.fill(HIST("mcTrack_Phi_PtvsEta"), mcTrack.pt(), mcTrack.eta());
                        }
                    }// if has_daughter

                }else{
                    // entered particles
                    if(mcTrack.has_daughters()){
                        auto daughters = mcTrack.daughters_as<aod::McParticles>();
                        int countmu_nomother = 0;
                        for(auto daughter : daughters){
                            if(abs(daughter.pdgCode())==13){
                                countmu_nomother++;
                            }
                        };

                        if(countmu_nomother == 2){
                            histos.fill(HIST("mcTrack_Phi_pt"), mcTrack.pt());
                            histos.fill(HIST("mcTrack_Phi_phi"), mcTrack.phi());
                            histos.fill(HIST("mcTrack_Phi_eta"), mcTrack.eta());
                            histos.fill(HIST("mcTrack_Phi_PtvsEta"), mcTrack.pt(), mcTrack.eta());
                        }
                    }// end if has daughters

                }// end if has mothers
            }// end if pdgCode == 13
        }// end mcTracks for
        std::cout << "finish calGenPhi" << std::endl;
    }; // end void traceDecay


    template <typename TMuon>
    bool singleMuoncut(TMuon const& muon,bool QA){

        if(QA){
            histos.fill(HIST("muon_BeforCut_trackType"), muon.trackType());
            histos.fill(HIST("muon_BeforCut_pt"), muon.pt());
            histos.fill(HIST("muon_BeforCut_eta"), muon.eta());
            histos.fill(HIST("muon_BeforCut_phi"), muon.phi());
            histos.fill(HIST("muon_BeforCut_pDCA"), muon.pDca());
            histos.fill(HIST("muon_BeforCut_Rabs"), muon.rAtAbsorberEnd());
            histos.fill(HIST("muon_BeforCut_chi2"), muon.chi2MatchMCHMFT());
        }
        if(!(muon.trackType()==0)){
            return false;
        };

        if(!(-3.6 <= muon.eta() && muon.eta() <=-2.5)){
            return false;
        };

        
        if(!(17.6 <= muon.rAtAbsorberEnd() && muon.rAtAbsorberEnd() <= 89.5)){
                return false;
        };

        if(!(0 <= muon.pDca() && muon.pDca() <= 594.0)){
                return false;
        };

        
            if(26.5 <= muon.rAtAbsorberEnd() &&muon.rAtAbsorberEnd() <= 89.5){
                if(!(0 <= muon.pDca() && muon.pDca() <= 324.0)){
                    return false;
                };
            };

            if(17.6 <= muon.rAtAbsorberEnd() &&muon.rAtAbsorberEnd() <= 26.5){
                if(!(0 <= muon.pDca() && muon.pDca() <= 594.0)){
                    return false;
                };
            };
        
        if(!(muon.chi2MatchMCHMFT() < fmatchingchi2cut)){
            return false;
        };
        if(QA){
            histos.fill(HIST("muon_AfterCut_trackType"), muon.trackType());
            histos.fill(HIST("muon_AfterCut_pt"), muon.pt());
            histos.fill(HIST("muon_AfterCut_eta"), muon.eta());
            histos.fill(HIST("muon_AfterCut_phi"), muon.phi());
            histos.fill(HIST("muon_AfterCut_pDCA"), muon.pDca());
            histos.fill(HIST("muon_AfterCut_Rabs"), muon.rAtAbsorberEnd());
            histos.fill(HIST("muon_AfterCut_chi2"), muon.chi2MatchMCHMFT());
        }
        return true;
    };

    template <typename TMuon,typename TMCTracks>
    bool mcMotherCut(TMuon const& muon,TMCTracks const& /*mcTracks*/,bool QA){
        if(!(muon.has_mcParticle())){
            return false;
        }else{
            auto mcMuon = muon.template mcParticle_as<TMCTracks>();
            auto mothers = mcMuon.template mothers_as<TMCTracks>();
            int countTarget = 0;
            for(auto mother : mothers){
                if(mother.pdgCode() == targetPdgCode){
                    if(QA){
                        histos.fill(HIST("muon_mother"),mother.pdgCode());
                    }
                    countTarget++;
                }
            }
            // mother <- phi only
            if(countTarget >= 1){
                return true;
            }else{
                return false;
            }
        }
    }

    enum MuonExtrapolation {
    // Index used to set different options for Muon propagation
    kToVertex = 0, // propagtion to vertex by default
    kToDCA,
    kToRabs
    };

    template <typename TEvents, typename TMuons,typename TMFTs,typename TMCTracks>
    void calRecoPhi(BCsWithTimestamps const& /*bcs*/,
                    TEvents const& events,
                    aod::FwdTrackAssoc const& muonAssocs,
                    TMuons const& /*muons*/,
                    TMFTs const& /*mfts*/,
                    TMCTracks const& mcTracks
                    ){
        o2::mch::TrackExtrap::setField();
        
        double refitpx = 0;
        double refitpy = 0;
        double refitpz = 0;
        double refitpt = 0;
        double refiteta = 0;
        double refitphi = 0;
        double refitsign = 0;
                
        double refitpx1 = 0;
        double refitpy1 = 0;
        double refitpz1 = 0;
        double refitpt1 = 0;
        double refiteta1 = 0;
        double refitphi1 = 0;
        double refitsign1 = 0;

        double refitpx2 = 0;
        double refitpy2 = 0;
        double refitpz2 = 0;
        double refitpt2 = 0;
        double refiteta2 = 0;
        double refitphi2 = 0;
        double refitsign2 = 0;

        for(auto assoc : muonAssocs){
            auto muon = assoc.template fwdtrack_as<TMuons>();
            for(auto event : events){
                if(muon.collisionId() == event.globalIndex()){
                    //propagate muon and global refit
                    if(muon.has_matchMCHTrack()&& muon.trackType()==0){
                        auto fwdbc = event.template bc_as<BCsWithTimestamps>();
                        initCCDB(fwdbc);
                        auto muontrack = muon.template matchMCHTrack_as<TMuons>();
                        auto mfttrack = muon.template matchMFTTrack_as<MyMFTTracks>();
                        //propagate muontrack to dca
                        o2::dataformats::GlobalFwdTrack propmuon = VarManager::PropagateMuon(muontrack, event, kToDCA);

                        //Refit (use MUON track p and MFT track eta,phi)
                        refitpx = propmuon.getP() * sin(M_PI / 2 - atan(mfttrack.tgl())) * cos(mfttrack.phi());
                        refitpy = propmuon.getP() * sin(M_PI / 2 - atan(mfttrack.tgl())) * sin(mfttrack.phi());
                        refitpz = propmuon.getP() * cos(M_PI / 2 - atan(mfttrack.tgl()));
                        
                            refitpt = std::sqrt(std::pow(refitpx, 2) + std::pow(refitpy, 2));
                            refiteta = mfttrack.eta();
                            refitphi =  mfttrack.phi();
                            refitsign = muon.sign();
                    }
                }
            }

            if(singleMuoncut(muon,true)){
                if(mcMotherCut(muon,mcTracks,true)){
                    // mother <- phi only
                    histos.fill(HIST("muon_pt"), refitpt);
                    histos.fill(HIST("muon_eta"), refiteta);
                    histos.fill(HIST("muon_phi"), refitphi);
                }
            };// end has_mcParticles
        }

        double MuonMass = 0.1056584;//GeV/c
        double dimuonMass = 0;
        //same event pairing
        for (auto event : events){
        //ToDo: need to use groupedAssocs , but I can't
             //auto groupedAssocs = muonAssocs.sliceBy(perCollision, event.globalIndex());
             //std::cout << "groupedAssocs.size() = " <<groupedAssocs.size() << std::endl;
            for (auto& [a1, a2] : o2::soa::combinations(muonAssocs, muonAssocs)) {
                if(!(a1.collisionId() == event.globalIndex() && a2.collisionId() == event.globalIndex())){
                    continue;
                }else{
                    auto t1 = a1.template fwdtrack_as<TMuons>();
                    auto t2 = a2.template fwdtrack_as<TMuons>();
                    if (t1.matchMCHTrackId() == t2.matchMCHTrackId()){
                        continue;
                    }
                    if (t1.matchMFTTrackId() == t2.matchMFTTrackId()){
                        continue;
                    }
                    if(!(singleMuoncut(t1,false))){
                        continue;
                    }
                    if(!(singleMuoncut(t2,false))){
                        continue;
                    }
                    if(!(mcMotherCut(t1,mcTracks,false))){
                        continue;
                    }
                    if(!(mcMotherCut(t2,mcTracks,false))){
                        continue;
                    }
                    auto mcTrack1 = t1.template mcParticle_as<TMCTracks>();
                    auto mcTrack2 = t2.template mcParticle_as<TMCTracks>();

                    auto mcMothers1 = mcTrack1.template mothers_as<TMCTracks>();
                    auto mcMothers2 = mcTrack2.template mothers_as<TMCTracks>();

                    bool sameMother = false;
                    for(auto mcMother1 : mcMothers1){
                        for(auto mcMother2 : mcMothers2){
                            if(mcMother1.pdgCode()==targetPdgCode && mcMother2.pdgCode() == targetPdgCode){
                                if(mcMother1.globalIndex() == mcMother2.globalIndex()){
                                    //std::cout << "smae global index: t1 t2 mother" << std::endl;
                                    sameMother = true;
                                }
                            }
                        }
                    }
                    if(!(sameMother)){
                        continue;
                    }
                    //calculate track1(refit)
                    //-------------------------------------
                    if(t1.has_matchMCHTrack() && t1.trackType()==0){
                        auto fwdbc1 = event.template bc_as<BCsWithTimestamps>();
                        initCCDB(fwdbc1);
                        auto muontrack1 = t1.template matchMCHTrack_as<TMuons>();
                        auto mfttrack1 = t1.template matchMFTTrack_as<MyMFTTracks>();
                        
                        //propagate muontrack to dca
                        o2::dataformats::GlobalFwdTrack propmuon1 = VarManager::PropagateMuon(muontrack1, event, kToDCA);

                        //Refit (use MUON track p and MFT track eta,phi)
                        refitpx1 = propmuon1.getP() * sin(M_PI / 2 - atan(mfttrack1.tgl())) * cos(mfttrack1.phi());
                        refitpy1 = propmuon1.getP() * sin(M_PI / 2 - atan(mfttrack1.tgl())) * sin(mfttrack1.phi());
                        refitpz1 = propmuon1.getP() * cos(M_PI / 2 - atan(mfttrack1.tgl()));
                        
                        refitpt1 = std::sqrt(std::pow(refitpx1, 2) + std::pow(refitpy1, 2));
                        refiteta1 = mfttrack1.eta();
                        refitphi1 =  mfttrack1.phi();
                        refitsign1 = t1.sign();
                    }
                    ////----------------------------------

                    //calculate track1(refit)
                    //-------------------------------------
                    if(t2.has_matchMCHTrack()&& t2.trackType()==0){
                        auto muontrack2 = t2.template matchMCHTrack_as<TMuons>();
                        auto mfttrack2 = t2.template matchMFTTrack_as<MyMFTTracks>();
                        //propagate muontrack to dca
                        o2::dataformats::GlobalFwdTrack propmuon2 = VarManager::PropagateMuon(muontrack2, event, kToDCA);

                        //Refit (use MUON track p and MFT track eta,phi)
                        refitpx2 = propmuon2.getP() * sin(M_PI / 2 - atan(mfttrack2.tgl())) * cos(mfttrack2.phi());
                        refitpy2 = propmuon2.getP() * sin(M_PI / 2 - atan(mfttrack2.tgl())) * sin(mfttrack2.phi());
                        refitpz2 = propmuon2.getP() * cos(M_PI / 2 - atan(mfttrack2.tgl()));
                        
                        refitpt2 = std::sqrt(std::pow(refitpx2, 2) + std::pow(refitpy2, 2));
                        refiteta2 = mfttrack2.eta();
                        refitphi2 =  mfttrack2.phi();
                        refitsign2 = t2.sign();
                    }
                    ////----------------------------------

                    //calculate dimuon invariant mass
                    ROOT::Math::PtEtaPhiMVector v1(refitpt1, refiteta1, refitphi1, MuonMass);
                    ROOT::Math::PtEtaPhiMVector v2(refitpt2, refiteta2, refitphi2, MuonMass);
                    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                    dimuonMass = v12.M();

                    
                    if(refitsign1 * refitsign2 > 0){
                        if(refitsign1 > 0){
                            histos.fill(HIST("dimuon_mass_SEPP"),dimuonMass);
                        }else if(refitsign1 < 0){
                            histos.fill(HIST("dimuon_mass_SEMM"),dimuonMass);
                        }
                    }else if(refitsign1 * refitsign2 < 0){
                        histos.fill(HIST("dimuon_mass_SEPM"),dimuonMass);
                        histos.fill(HIST("dimuon_Phi_PtvsEta"),v12.pt(),v12.eta());
                    }
                }
                std::cout << "++++++++++++++++  combination end  ++++++++++++++++++++++++++" << std::endl;
                //VarManager::FillPairVertexing(event, t1, t2, true);
            }//end combination
        }//end event

        std::cout << "finish calRecoPhi" << std::endl;
    };

    // main process
    void process(
                BCsWithTimestamps const& bcs,
                MyEvents const& events,
                aod::FwdTrackAssoc const& muonAssocs,
                MyMuonsWithCov const& muons,
                MyMFTTracks const& mfts,
                aod::McCollisions const& mcEvents,
                aod::McParticles const& mcTracks){
        calGenPhi(mcEvents,mcTracks);
        calRecoPhi(bcs,events,muonAssocs,muons,mfts,mcTracks);
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& config) {
  return WorkflowSpec{
      adaptAnalysisTask<DimuonAccEff>(config),
  };
}
