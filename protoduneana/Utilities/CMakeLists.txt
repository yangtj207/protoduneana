cet_find_library( DUNEARTSUPPORT NAMES dune_ArtSupport PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
cet_find_library( PROTODUNEDATAUTILS NAMES ProtoDUNEDataUtils PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
cet_find_library( DUNEOBJLD NAMES dune_DuneObj PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )

art_make( LIBRARY_NAME ProtoDUNEUtilities
	  LIB_LIBRARIES
                        DUNEARTSUPPORT
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larreco_RecoAlg
                        larreco_Calorimetry
                        larreco_RecoAlg_PMAlg
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
                        lardata_ArtDataHelper
                        lardata_Utilities
                        larsim_Simulation 
                        lardataobj_Simulation
                        PROTODUNEDATAUTILS
                        geant4reweight::ReweightBaseLib
                        geant4reweight::PropBaseLib
                        DUNEOBJLD
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        nug4_ParticleNavigation 
                        nusimdata_SimulationBase
                        ${CLHEP}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
	  MODULE_LIBRARIES
          ProtoDUNEUtilities
)

install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(FitUtils)
