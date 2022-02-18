
art_make(BASENAME_ONLY
  MODULE_LIBRARIES
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common
  canvas
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  cetlib 
  cetlib_except
  lardataobj::RecoBase
  lardata::Utilities
  larcorealg::Geometry
  lardataobj::RawData
  nusimdata::SimulationBase
  larsim::MCCheater_PhotonBackTrackerService_service           
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  ROOT_BASIC_LIB_LIST
  )
if(DEFINED ENV{HEP_HPC_DIR})
  add_subdirectory(hdf)
endif()

install_headers()
install_fhicl()
install_source()
