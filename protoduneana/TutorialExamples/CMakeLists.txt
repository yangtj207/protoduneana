cet_build_plugin(ClusterCounter art::module
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter2 art::module
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter3 art::module
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter4 art::module
                        lardataalg_DetectorInfo
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ClusterCounter5 art::module
                        lardataalg_DetectorInfo
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ReadSpacePointAndCnn art::module
			lardata_ArtDataHelper
			lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(ReadBeamInfo art::module
      nusimdata_SimulationBase
      larcorealg_Geometry
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
      ${ART_ROOT_IO_TFILESERVICE_SERVICE}
      ${MF_MESSAGELOGGER}
      ROOT_BASIC_LIB_LIST
      BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
