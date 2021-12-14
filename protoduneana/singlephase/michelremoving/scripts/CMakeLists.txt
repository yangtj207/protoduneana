cet_make_exec( calculate_calibration_constants
               SOURCE calconst.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               )

cet_make_exec( fitlangaumpv
               SOURCE fitlangaumpv.C
               LIBRARIES
               ROOT_BASIC_LIB_LIST
               )

install_scripts()
