list(APPEND TrueLibs dunecore_ArtSupport ProtoDUNEDataUtils dunecore_DuneObj)
list(APPEND AliasLibs DUNEARTSUPPORT PROTODUNEDATAUTILS DUNEOBJLD)

foreach(lib IN ZIP_LISTS TrueLibs AliasLibs)
  if(DEFINED ENV{DUNETPC_FQ_DIR})
    message("DUNETPC_FQ_DIR defined. Using this to search for library ${lib_0}")
    cet_find_library( ${lib_1} NAMES ${lib_0} PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
  else()
    message("DUNETPC_FQ_DIR undefined. Expect we are building next to dunetpc. Aliasing library ${lib_0} to ${lib_1}")
    #add_library(${lib_1} ALIAS ${lib_0})
    set(${lib_1} ${lib_0})
  endif()
endforeach()

#cet_find_library( DUNEARTSUPPORT NAMES dunecore_ArtSupport PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
#cet_find_library( PROTODUNEDATAUTILS NAMES ProtoDUNEDataUtils PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
#cet_find_library( DUNEOBJLD NAMES dunecore_DuneObj PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )

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
                        duneprototypes::ProtoDUNEDataUtils
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
