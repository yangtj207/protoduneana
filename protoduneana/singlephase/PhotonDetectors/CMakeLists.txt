
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

cet_build_plugin(PDSPmatch art::module
                        lardataobj::RawData
                        lardataobj::RecoBase
                        larcorealg::Geometry
                        dune_raw_data::dune-raw-data_Overlays
                        artdaq_core::artdaq-core_Data
                        dunecore::Geometry
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib 
                        cetlib_except::cetlib_except
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
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support
                        ROOT::Core
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        canvas::canvas
                        art::Persistency_Provenance
                        art::Utilities
                        messagefacility::MF_MessageLogger
                        cetlib::cetlib 
                        cetlib_except::cetlib_except
                        ROOT_BASIC_LIB_LIST
                        RCEDAMLIB
                        BASENAME_ONLY
)


install_fhicl()
install_source()
install_scripts()
