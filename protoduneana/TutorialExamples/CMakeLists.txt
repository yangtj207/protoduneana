cet_build_plugin(ClusterCounter art::module
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter2 art::module
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter3 art::module
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter4 art::module
                        lardataalg::DetectorInfo
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter5 art::module
                        lardataalg::DetectorInfo
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
			nusimdata::SimulationBase
			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ReadSpacePointAndCnn art::module
			lardata::ArtDataHelper
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			nusimdata::SimulationBase
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			messagefacility::MF_MessageLogger
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ReadBeamInfo art::module
      nusimdata::SimulationBase
      larcorealg::Geometry
      ART_FRAMEWORK_SERVICES_REGISTRY
      ART_ROOT_IO_TFILE_SUPPORT
      ROOT::Core
      ART_ROOT_IO_TFILESERVICE_SERVICE
      messagefacility::MF_MessageLogger
      ROOT_BASIC_LIB_LIST
      BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
