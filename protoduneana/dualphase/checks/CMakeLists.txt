cet_find_library( PROTODUNEDATAUTILS NAMES ProtoDUNEDataUtils PATHS ENV DUNETPC_FQ_DIR PATH_SUFFIXES lib NO_DEFAULT_PATH )

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
