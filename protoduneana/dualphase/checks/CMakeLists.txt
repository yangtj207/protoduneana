#cet_find_library( PROTODUNEDATAUTILS NAMES ProtoDUNEDataUtils PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
list(APPEND TrueLibs ProtoDUNEDataUtils)
list(APPEND AliasLibs PROTODUNEDATAUTILS)

foreach(lib IN ZIP_LISTS TrueLibs AliasLibs)
  if(DEFINED ENV{DUNETPC_FQ_DIR})
    message("DUNETPC_FQ_DIR defined. Using this to search for library ${lib_0}")
    cet_find_library( ${lib_1} NAMES ${lib_0} PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )
  else()
    message("DUNETPC_FQ_DIR undefined. Expect we are building next to dunetpc. Aliasing library ${lib_0} to ${lib_1}")
    set(${lib_1} ${lib_0})
  endif()
endforeach()

cet_build_plugin(		HitNormCheck art::module
			lardataobj_RecoBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
                        BASENAME_ONLY
)

cet_build_plugin(		CosmicsdQdx art::module
			lardataobj_RecoBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
			PROTODUNEDATAUTILS
			ProtoDUNEUtilities
                        BASENAME_ONLY
)

cet_build_plugin(		TrackHitInfo art::module
			lardataobj_RecoBase
                        lardata_ArtDataHelper
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
			PROTODUNEDATAUTILS
			ProtoDUNEUtilities
                        BASENAME_ONLY
)

install_fhicl()
install_source()
