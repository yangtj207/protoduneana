cet_make_exec( dEdX_calibration 
               SOURCE protoDUNE_dEdx_calib.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               CaloTools
	       cetlib 
	       cetlib_except
	       fhiclcpp
               messagefacility::MF_MessageLogger
               )
