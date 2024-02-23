cet_make_exec( plotdeltaz
               SOURCE plotdeltaz.C
               LIBRARIES
               ROOT::Core ROOT::Tree ROOT::Hist
               fhiclcpp::fhiclcpp
               )

install_scripts()
install_headers()
install_fhicl()
