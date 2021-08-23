cet_make_exec( make_x_correction
               SOURCE protoDUNE_X_calib.C
               LIBRARIES
               ${ROOT_BASIC_LIB_LIST}
               )

cet_make_exec( validate_calibration
               SOURCE protoDUNE_validate_calib.C
               LIBRARIES
               ${ROOT_BASIC_LIB_LIST}
               )
