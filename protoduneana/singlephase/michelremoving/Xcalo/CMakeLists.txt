find_library(CaloTools ../scripts) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - prefer find_package() with custom Find<pkg>.cmake for transitivity, relocatability
cet_make_exec( make_x_correction
               SOURCE protoDUNE_X_calib.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               CaloTools              
               )
install_headers()
install_fhicl()
install_source()
install_scripts()

