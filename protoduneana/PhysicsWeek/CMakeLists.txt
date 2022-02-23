cet_build_plugin(MichelHitCounter art::module
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              nusimdata::SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
	      ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

cet_build_plugin(VtxScore art::module
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              nusimdata::SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
	      ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)


cet_build_plugin(MichelEventSelection art::module
              larcorealg::Geometry
              larreco::Calorimetry
              larreco::RecoAlg_PMAlg
              larsim::MCCheater_BackTrackerService_service
              larsim::MCCheater_ParticleInventoryService_service
              lardataobj::RecoBase
              lardata::ArtDataHelper
              nusimdata::SimulationBase
              ART_FRAMEWORK_CORE
              ART_FRAMEWORK_PRINCIPAL
              ART_FRAMEWORK_SERVICES_REGISTRY
              ART_ROOT_IO_TFILE_SUPPORT
	      ROOT::Core
              ART_ROOT_IO_TFILESERVICE_SERVICE
              art::Persistency_Provenance
              messagefacility::MF_MessageLogger
              cetlib::cetlib
              ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
