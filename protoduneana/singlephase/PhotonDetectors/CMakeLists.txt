
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

cet_build_plugin(PDSPmatch art::module
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larcorealg::Geometry
                        dune_raw_data::dune-raw-data_Overlays
                        artdaq_core::artdaq-core_Data
                        dunecore::Geometry
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT 
			ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Provenance
                        messagefacility::MF_MessageLogger
                        cetlib
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)

cet_build_plugin(PDSPmatchmc art::module
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larcorealg::Geometry
                        dune_raw_data::dune-raw-data_Overlays
                        artdaq_core::artdaq-core_Data
                        dunecore::Geometry
                        nusimdata::SimulationBase
                        larsim::MCCheater_PhotonBackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			ART_PERSISTENCY_PROVENANCE
                        messagefacility::MF_MessageLogger
                        cetlib
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)


install_fhicl()
install_source()
install_scripts()
