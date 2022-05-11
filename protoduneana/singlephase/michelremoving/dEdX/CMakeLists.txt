cet_make_exec( dEdX_calibration 
               SOURCE protoDUNE_dEdx_calib.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               CaloTools
	       cetlib::cetlib 
	       cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               messagefacility::MF_MessageLogger
               )

art_make_library( LIBRARY_NAME protodune_dedx_fit
		  SOURCE ProtoDUNEdEdXFitter.cxx 
		  LIBRARIES ${ART_FRAMEWORK_CORE}
	       		    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_ROOT_IO_TFILE_SUPPORT} 
			    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			    ${ART_PERSISTENCY_COMMON}
			    canvas
			    art_Persistency_Provenance
			    cetlib cetlib_except
			    ${ART_UTILITIES}
			    ${MF_MESSAGELOGGER}
			    ${FHICLCPP}
			    ${CETLIB}
			    ROOT_BASIC_LIB_LIST
                            CaloTools
                            
			    ${ROOT_ROOFITCORE}
			    ${ROOT_ROOFIT}
			    ${ROOT_ROOSTATS}
			    ${ROOT_HISTFACTORY}
			    ${ROOT_FOAM}
			    ${ROOT_GEOM}
			    ${ROOT_MINUIT}
			    ${ROOT_MINUIT2}
			    ${ROOT_CORE}
			    ${ROOT_MATHCORE}
	    )

art_make_exec( NAME runProtoDUNEdEdXFit
	       SOURCE runProtoDUNEdEdXFit.c++
	       LIBRARIES protodune_dedx_fit 
	    )
install_scripts()
