install_fhicl()
install_scripts()

INSTALL(FILES asimov.cfg
        DESTINATION cfg_files) ### MIGRATE-CHANGED (migrate-3.22.02) - ${product}/+${version}/* -> ""
