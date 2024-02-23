art_make(BASENAME_ONLY
  MODULE_LIBRARIES
  cetlib cetlib_except
  larsim::Simulation
  lardataobj::Simulation
  lardataobj::RecoBase
  lardata::Utilities
  lardata::ArtDataHelper
  larcorealg::Geometry
  nusimdata::SimulationBase
  larsim::MCCheater_PhotonBackTrackerService_service           
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  ROOT::Core ROOT::Tree ROOT::Hist
  ProtoDUNEUtilities
  geant4reweight::ReweightBaseLib
  geant4reweight::PropBaseLib
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  canvas::canvas
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  hep_hpc_hdf5
  hep_hpc_Utilities
  ${HDF5_LIBRARIES}
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types
install_headers()
install_fhicl()
install_source()
