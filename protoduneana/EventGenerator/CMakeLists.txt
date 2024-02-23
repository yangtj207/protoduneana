
include_directories( $ENV{IFDHC_FQ_DIR}/inc ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible
art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
	  dunesim::EventGenerator_ProtoDUNEbeamDataProducts
          larcoreobj::SummaryData
          larcorealg::Geometry
          larcore::Geometry_Geometry_service
          nurandom::RandomUtils_NuRandomService_service
          larsim::Simulation lardataobj_Simulation
	  lardataobj::RecoBase
          nusimdata_SimulationBase
          art::Framework_Core
          art::Framework_Principal
          art::Framework_Services_Registry
          art::Framework_Services_Optional_RandomNumberGenerator_service
          art::Persistency_Common
          art::Persistency_Provenance
          art::Utilities 
	  canvas::canvas
          messagefacility::MF_MessageLogger
          hep_concurrency
          fhiclcpp::fhiclcpp
          cetlib::cetlib 
	  cetlib_except::cetlib_except
          ${CLHEP}
          ROOT::Core ROOT::Tree ROOT::Hist
          ${SQLITE3}
          ${ROOT_EG}
          ifdh::ifdh 
	  duneprototypes::ProtoDUNEDataUtils
          ProtoDUNEUtilities
          art_root_io::tfile_support
	  ROOT::Core
          art_root_io::TFileService_service
 ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
