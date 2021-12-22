cet_make_exec( calculate_calibration_constants
               SOURCE calconst.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               )

cet_make_library( LIBRARY_NAME CaloTools
                  SOURCE
                  LanGausFit.C
                  LIBRARIES
                  ROOT_BASIC_LIB_LIST
                  )

cet_make_exec( fitlangaumpv
               SOURCE fitlangaumpv.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               CaloTools
               )

cet_make_exec( murecombfit
               SOURCE murecombfit.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               ROOT_LIBRARIES
	       cetlib cetlib_except
	       FHICLCPP
               MF_MESSAGELOGGER
               CaloTools
               )

install_scripts()
install_headers()
