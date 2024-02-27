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

cmake_minimum_required (VERSION 3.20 FATAL_ERROR) ### MIGRATE-CHANGED (migrate-3.22.02) - VERSION 3.19 -> 3.20
# cmake_policy(VERSION 3.19) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.22.02) - Uncomment to preserve compatibility with older CMake versions

find_package(cetmodules REQUIRED) ### MIGRATE-ADDED (migrate-3.22.02) - find_package(cetmodules) must precede project()
project(protoduneana LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.83.01d00)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

# Find our own CMake modules, and tell everyone else where to find them.
cet_cmake_module_directories(Modules BINARY)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined
)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - add args: REPORT_THRESHOLD VERBOSE

# these are minimum required versions, not the actual product versions
find_package( Boost REQUIRED ) # Needs to be first to cover for a missing transitive ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
                  # dependency from artdaq_core.
find_package( Geant4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package(ROOT REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)

find_package( dunecore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunecore ...) -> find_package(dunecore REQUIRED ...)
find_package( duneprototypes REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(duneprototypes ...) -> find_package(duneprototypes REQUIRED ...)
find_package( dunecalib REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunecalib ...) -> find_package(dunecalib REQUIRED ...)
find_package( duneopdet REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(duneopdet ...) -> find_package(duneopdet REQUIRED ...)
find_package( dunesim REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunesim ...) -> find_package(dunesim REQUIRED ...)
find_package( larcoreobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcoreobj ...) -> find_package(larcoreobj REQUIRED ...)
find_package( larcore  REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcore ...) -> find_package(larcore REQUIRED ...)
find_package( larevt REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larevt ...) -> find_package(larevt REQUIRED ...)
find_package( lardataobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataobj ...) -> find_package(lardataobj REQUIRED ...)
find_package( lardataalg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataalg ...) -> find_package(lardataalg REQUIRED ...)
find_package( lardata  REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardata ...) -> find_package(lardata REQUIRED ...)
find_package( larsim REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larsim ...) -> find_package(larsim REQUIRED ...)
find_package( larpandora  REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larpandora ...) -> find_package(larpandora REQUIRED ...)
#find_package( pandora REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(pandora ...) -> find_package(pandora REQUIRED ...)
find_package( larreco REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larreco ...) -> find_package(larreco REQUIRED ...)
find_package( larrecodnn REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larrecodnn ...) -> find_package(larrecodnn REQUIRED ...)
find_package( larana REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larana ...) -> find_package(larana REQUIRED ...)
find_package( nutools REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nutools ...) -> find_package(nutools REQUIRED ...)
find_package( nug4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nug4 ...) -> find_package(nug4 REQUIRED ...)
find_package( nugen REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nugen ...) -> find_package(nugen REQUIRED ...)
find_package( GENIE REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(genie ...) -> find_package(genie REQUIRED ...)
find_package( nurandom REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nurandom ...) -> find_package(nurandom REQUIRED ...)
find_package( art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art ...) -> find_package(art REQUIRED ...)
find_package( art_root_io REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art_root_io ...) -> find_package(art_root_io REQUIRED ...)
find_package( CLHEP REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(clhep ...) -> find_package(clhep REQUIRED ...)
find_package( messagefacility REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(messagefacility ...) -> find_package(messagefacility REQUIRED ...)
find_package( PostgreSQL REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(postgresql ...) -> find_package(postgresql REQUIRED ...)
find_package( dunepdlegacy REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunepdlegacy ...) -> find_package(dunepdlegacy REQUIRED ...)
find_package( artdaq_core REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(artdaq_core ...) -> find_package(artdaq_core REQUIRED ...)
find_package( TRACE REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(TRACE ...) -> find_package(TRACE REQUIRED ...)
find_package( ifdh_art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(ifdh_art ...) -> find_package(ifdh_art REQUIRED ...)
find_package( ifdhc REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(ifdhc ...) -> find_package(ifdhc REQUIRED ...)
find_package( systematicstools REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(systematicstools ...) -> find_package(systematicstools REQUIRED ...)
find_package( nusystematics REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(nusystematics ...) -> find_package(nusystematics REQUIRED ...)
find_package( geant4reweight REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(geant4reweight ...) -> find_package(geant4reweight REQUIRED ...)
#find_package( valgrind REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(valgrind ...) -> find_package(valgrind REQUIRED ...)

if(DEFINED ENV{HEP_HPC_DIR} )
  find_package( hdf5 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(hdf5 ...) -> find_package(hdf5 REQUIRED ...)
  find_package( hep_hpc REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(hep_hpc ...) -> find_package(hep_hpc REQUIRED ...)
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

# packaging utility
cet_cmake_config() ### MIGRATE-ADDED (migrate-3.22.02) - required cmd
