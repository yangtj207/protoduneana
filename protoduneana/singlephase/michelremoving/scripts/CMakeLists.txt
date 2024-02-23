cet_make_exec( validate_calibration
               SOURCE protoDUNE_validate_calib.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
	       ROOT::Physics
               fhiclcpp::fhiclcpp
               )

cet_make_exec( calculate_calibration_constants
               SOURCE calconst.C
               LIBRARIES
	       ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
               )

cet_make_library( LIBRARY_NAME CaloTools
                  SOURCE
                  LanGausFit.C CaloUtils.cxx
                  LIBRARIES
                  ROOT::Core ROOT::Tree ROOT::Hist
                  )

cet_make_exec( fitlangaumpv
               SOURCE fitlangaumpv.C
               LIBRARIES
	       ROOT::Core ROOT::Tree ROOT::Hist ROOT::Gpad
               CaloTools
               )

cet_make_exec( murecombfit
               SOURCE murecombfit.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
	       ROOT::Minuit
	       cetlib::cetlib cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               CaloTools
               )

cet_make_exec( plotmurecombfitresult
               SOURCE plotmurecombfitresult.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
               ROOT::Minuit
               cetlib::cetlib cetlib_except::cetlib_except
               fhiclcpp::fhiclcpp
               CaloTools
               )

cet_make_exec( resrangefit
               SOURCE resrangefit.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               ROOT_LIBRARIES
	       ROOT::Minuit
               cetlib::cetlib cetlib_except::cetlib_except
	       fhiclcpp::fhiclcpp
               CaloTools
               )

cet_make_exec( ConvertXCorrection
               SOURCE ConvertXCorrection.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               )

cet_make_exec( ConvertYZCorrection
               SOURCE ConvertYZCorrection.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               )

install_scripts()
install_headers()
