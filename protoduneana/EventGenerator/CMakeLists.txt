
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
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
          ART_FRAMEWORK_CORE
          ART_FRAMEWORK_PRINCIPAL
          ART_FRAMEWORK_SERVICES_REGISTRY
          ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
          art::Persistency_Common
          art::Persistency_Provenance
          art::Utilities 
	  canvas
          messagefacility::MF_MessageLogger
          hep_concurrency
          fhiclcpp
          cetlib 
	  cetlib_except
          ${CLHEP}
          ROOT_BASIC_LIB_LIST
          ${SQLITE3}
          ${ROOT_EG}
          ifdh::ifdh 
	  duneprototypes::ProtoDUNEDataUtils
          ProtoDUNEUtilities
          ART_ROOT_IO_TFILE_SUPPORT
	  ROOT::Core
          ART_ROOT_IO_TFILESERVICE_SERVICE
 )


install_headers()
install_fhicl()
install_source()
