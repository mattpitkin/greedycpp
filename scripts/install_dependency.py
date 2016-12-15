#!/usr/bin/env python
###############################################################################
# Author: Ka Wa TSANG
# Data  : 20161202
# Description :
#   A script to download and install packages (tar.gz).
#   It creates a for loop over the variable "Packages", then do the following.
#   1) wget "link" -P "INSTALL_SRC"
#   2) cd "INSTALL_SRC"
#   3) tar zxvf "NAME.tar.gz"
#   4) cd "INSTALL_SRC/NAME.tar.gz"
#   5) ./configure --prefix="INSTALL_DIR" "OPTIONS" && make && make install
###############################################################################

import os 

# User input (Edited it if needed) 
INSTALL_SRC = "${HOME}/install/src"   # Path to save the installation file
INSTALL_DIR = "${HOME}/install/opt"   # Path for configuring prefix

os.system( "mkdir -p " + INSTALL_SRC )
os.system( "mkdir -p " + INSTALL_DIR )

###############################################################################
# Package information
# [ name (w/o tar.gz) , link , configure option ]
cmake     = [ "cmake-3.7.0"        , "https://cmake.org/files/v3.7/cmake-3.7.0.tar.gz"                            , ""                 ]
libconfig = [ "libconfig-1.5"      , "http://www.hyperrealm.com/libconfig/libconfig-1.5.tar.gz"                   , ""                 ]
gsl       = [ "gsl-2.2.1"          , "http://nl.mirror.babylon.network/gnu/gsl/gsl-2.2.1.tar.gz"                  , ""                 ]
openmpi   = [ "openmpi-2.0.1"      , "https://www.open-mpi.org/software/ompi/v2.0/downloads/openmpi-2.0.1.tar.gz" , "--enable-mpi-cxx" ]
hdf5      = [ "hdf5-1.10.0-patch1" , "ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.10.0-patch1.tar.gz"          , ""                 ]

# TODO: Add installation steps for git repository.
cnpy      = [ "https://github.com/rogersce/cnpy.git" ]

# Package you want to install
# Packages = [ cmake , libconfig , gsl , openmpi , hdf5 ]
Packages = [ cmake , libconfig , gsl , openmpi , hdf5 ]

###############################################################################

# Main
for pack in Packages:
  # Download packages if needed
  print " --- "
  print "Checking %s.tar.gz ..." % (pack[0])
  if os.path.exists( "%s/%s.tar.gz" % (INSTALL_SRC,pack[0]) ):
    print "%s.tar.gz exists." % (pack[0])
  else:
    print "%s.tar.gz doesnt exist. Downloading ..." % (pack[0])
    os.system( "wget %s -P %s" % (pack[1],INSTALL_SRC) )
    print "%s.tar.gz is downloaded" % (pack[0])
  
  os.chdir( INSTALL_SRC )
  print "Untarring %s.tar.gz ..." % (pack[0])
  os.system( "tar zxvf %s.tar.gz" % (pack[0]) )

  os.chdir( "%s/%s" % (INSTALL_SRC,pack[0]) )
  print "Configuring and making %s" % (pack[0])
  os.system( "./configure --prefix=%s %s && make && make install" % (INSTALL_DIR,pack[2]) )
 