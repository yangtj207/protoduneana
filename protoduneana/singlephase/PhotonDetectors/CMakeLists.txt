cet_find_library( DUNEGEOMETRY NAMES dune_Geometry PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )

cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

cet_build_plugin(PDSPmatch art::module
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        dune_raw_data::dune-raw-data_Overlays
                        artdaq_core::artdaq-core_Data
                        DUNEGEOMETRY
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)

cet_build_plugin(PDSPmatchmc art::module
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        dune_raw_data::dune-raw-data_Overlays
                        artdaq_core::artdaq-core_Data
                        DUNEGEOMETRY
                        nusimdata_SimulationBase
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_PROVENANCE}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)


install_fhicl()
install_source()
install_scripts()
