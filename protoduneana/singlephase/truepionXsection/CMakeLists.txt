
art_make(BASENAME_ONLY
  MODULE_LIBRARIES
  lardataobj::RecoBase
  lardata::Utilities
  larcorealg::Geometry
  nusimdata::SimulationBase
  larsim::MCCheater_PhotonBackTrackerService_service           
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  ROOT_BASIC_LIB_LIST
  duneprototypes::ProtoDUNEDataUtils
  ProtoDUNEUtilities
  larsim::Simulation
  lardataobj::Simulation
  lardataobj::RawData
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  art_root_io
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities 
  canvas
  messagefacility::MF_MessageLogger
  cetlib 
  cetlib_except
  )

cet_remove_compiler_flags(CXX -pedantic)

install_headers()
install_fhicl()
install_source()
