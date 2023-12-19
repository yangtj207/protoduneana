art_make(BASENAME_ONLY
  TOOL_LIBRARIES
  art::Framework_Core
  larcorealg::Geometry
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  ROOT_BASIC_LIB_LIST     
  hep_hpc_hdf5
  ${HDF5_LIBRARIES}

  )

install_headers()
install_fhicl()
install_source()
