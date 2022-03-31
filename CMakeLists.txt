# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================

cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(protoduneana VERSION 09.47.00.00 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

# Find our own CMake modules, and tell everyone else where to find them.
cet_cmake_module_directories(Modules BINARY)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_boost( ) # Needs to be first to cover for a missing transitive
                  # dependency from artdaq_core.
find_ups_geant4( )
find_ups_root()

find_ups_product( dunecore )
find_ups_product( duneprototypes )
find_ups_product( dunecalib )
find_ups_product( duneopdet )
find_ups_product( dunesim )
find_ups_product( larcoreobj )
find_ups_product( larcore  )
find_ups_product( larevt )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata  )
find_ups_product( larsim )
find_ups_product( larpandora  )
find_ups_product( pandora )
find_ups_product( larreco )
find_ups_product( larana )
find_ups_product( nutools )
find_ups_product( nug4 )
find_ups_product( nugen )
find_ups_product( genie )
find_ups_product( nurandom )
find_ups_product( art )
find_ups_product( art_root_io)
find_ups_product( clhep )
find_ups_product( messagefacility )
find_ups_product( cetbuildtools )
find_ups_product( postgresql )
find_ups_product( dune_raw_data )
find_ups_product( dunepdsprce )
find_ups_product( artdaq_core )
find_ups_product( ifdh_art )
find_ups_product( ifdhc )
find_ups_product( systematicstools )
find_ups_product( nusystematics )
find_ups_product( geant4reweight )
if(DEFINED ENV{HEP_HPC_DIR} )
  find_ups_product( hdf5 )
  find_ups_product( hep_hpc )
endif()

if(DEFINED ENV{HDF5_DIR} )
  string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
  if( ${BTYPE_UC} MATCHES "DEBUG" )
    cet_find_library( HDF5 NAMES hdf5_debug PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
  else()
    cet_find_library( HDF5 NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
  endif()
endif()

# macros for dictionary and cet_build_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(protoduneana)
add_subdirectory(fcl)

# Unit tests.
#add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
