# Makefile for the CKDMIP software
#
# This makefile builds two executables: lblrtm2nc, which converts
# LBLRTM binary files to NetCDF4/HDF5, and ckdmip_lw, which performs
# longwave radiative transfer.  You may need to edit the configuration
# entries below for your compiler.

# --- CONFIGURATION ---

# First try using the nf-config utility to set the NETCDF_INCLUDE and
# NETCDF_LIB flags
HAVE_NFCONFIG := $(shell nf-config --version 2> /dev/null)

ifdef HAVE_NFCONFIG

$(info *** Using nf-config to obtain NetCDF flags)
NETCDF_INCLUDE = $(shell nf-config --fflags)
NETCDF_LIB     = $(shell nf-config --flibs)

else

$(info *** nf-config not found)
# If the Fortran NetCDF library has not been compiled into a system
# directory, you need to define NETCDF_INCLUDE and NETCDF_LIB to
# indicate the location of the NetCDF4 Fortran module (netcdf.mod) and
# the NetCDF4 static or shared library.  At ECMWF these are defined
# with "module load netcdf4" leading to them containing the following:
#
#NETCDF_INCLUDE = -I/usr/local/apps/netcdf4/4.6.2/GNU/7.3.0/include
#
#NETCDF_LIB = -L/usr/local/apps/netcdf4/4.6.2/GNU/7.3.0/lib -lnetcdff -lnetcdf_c++ -lnetcdf -Wl,-rpath,/usr/local/apps/netcdf4/4.6.2/GNU/7.3.0/lib

endif

# Fortran compiler executable
FC = gfortran

# Tell the compiler to accept C-preprocessor directives, and compile
# easy_netcdf.F90 with NetCDF4/HDF5 support
CPPFLAGS = -cpp -DNC_NETCDF4

# Basic flags such as where to write module files
BASICFLAGS =

# OpenMP flag
OMPFLAG = -fopenmp

# Warning flags: all except those that warn about unused stuff
WARNFLAGS = -Wall -Wno-unused-label -Wno-unused-dummy-argument -Wno-unused-variable

# Debugging flags, such as "-g" to store debugging symbols, and
# instructions to abort if certain floating-point exceptions occur
DEBUGFLAGS = -g -ffpe-trap=invalid,zero,overflow
OPTFLAGS = -O2 -march=native

# For debugging use the following arguments for gfortran
DEBUGFLAGS = -g -ffpe-trap=invalid,zero,overflow -fcheck=bounds -finit-real=snan
#OPTFLAGS = -O0

# Uncomment the following for PGI Fortran
#FC = pgf90
#CPPFLAGS = -Mpreprocess
#BASICFLAGS = -module mod 
#DEBUGFLAGS = -g
#WARNFLAGS = 
#OMPFLAG = -mp=nonuma
#OPTFLAGS = -O0


# --- CONSOLIDATE VARIABLES AND COMPILER ARGUMENTS ---

ifndef NETCDF_INCLUDE
$(warning Warning: NETCDF_INCLUDE is not defined - you may need to edit src/Makefile)
endif
ifndef NETCDF_LIB
$(warning Warning: NETCDF_LIB is not defined     - you may need to edit src/Makefile)
endif

export FCFLAGS = $(WARNFLAGS) $(BASICFLAGS) $(CPPFLAGS) \
	$(OPTFLAGS) $(DEBUGFLAGS) $(NETCDF_INCLUDE) $(OMPFLAG)
export LIBS = $(LDFLAGS) $(NETCDF_LIB) $(OMPFLAG)

BASESOURCES = easy_netcdf.F90 radiation_io.F90 yomlun_ifsaux.F90 \
	parkind1.F90 abor1.F90 interpolation.F90

RADLW_SOURCES = \
	ckdmip_lw.F90 longwave_config.F90 \
	longwave_fluxes.F90 planck_function.F90

RADSW_SOURCES = \
	ckdmip_sw.F90 shortwave_config.F90 \
	shortwave_fluxes.F90 
#shortwave_tools.F90

TOOL_SOURCES = ckdmip_tool.F90

CONVERT_SOURCES = ckdmip_convert.F90

LBL2NC_SOURCES = lblrtm2nc.F90

CONTINUUM_SOURCES = change_continuum.F90 caviar_continuum.F90 continuum_lut.F90

BASEOBJECTS := $(BASESOURCES:.F90=.o)

RADLW_OBJECTS   := $(BASEOBJECTS) $(RADLW_SOURCES:.F90=.o)
RADSW_OBJECTS   := $(BASEOBJECTS) $(RADSW_SOURCES:.F90=.o)
LBL2NC_OBJECTS  := $(BASEOBJECTS) $(LBL2NC_SOURCES:.F90=.o)
TOOL_OBJECTS    := $(BASEOBJECTS) $(TOOL_SOURCES:.F90=.o)
CONVERT_OBJECTS := $(BASEOBJECTS) $(CONVERT_SOURCES:.F90=.o)
CONTINUUM_OBJECTS:=$(BASEOBJECTS) $(CONTINUUM_SOURCES:.F90=.o)
OBJECTS := $(BASEOBJECTS) $(RADLW_SOURCES:.F90=.o) \
	$(RADSW_SOURCES:.F90=.o) $(LBL2NC_SOURCES:.F90=.o) \
	$(TOOL_SOURCES:.F90=.o) $(CONVERT_SOURCES:.F90=.o) \
	$(CONTINUUM_SOURCES:.F90=.o)

RADLW = ../bin/ckdmip_lw
RADSW = ../bin/ckdmip_sw
TOOL = ../bin/ckdmip_tool
LBL2NC = ../bin/lblrtm2nc
CONVERT = ../bin/ckdmip_convert
CONTINUUM = ../bin/change_continuum

PROGRAMS = $(RADLW) $(RADSW) $(LBL2NC) $(TOOL) $(CONVERT) $(CONTINUUM)

# --- BUILD TARGETS ---

all: $(PROGRAMS)

$(RADLW): $(RADLW_OBJECTS)
	$(FC) $(FCFLAGS) $(RADLW_OBJECTS) $(LIBS) -o $(RADLW)

$(RADSW): $(RADSW_OBJECTS)
	$(FC) $(FCFLAGS) $(RADSW_OBJECTS) $(LIBS) -o $(RADSW)

$(TOOL): $(TOOL_OBJECTS)
	$(FC) $(FCFLAGS) $(TOOL_OBJECTS) $(LIBS) -o $(TOOL)

$(CONVERT): $(CONVERT_OBJECTS)
	$(FC) $(FCFLAGS) $(CONVERT_OBJECTS) $(LIBS) -o $(CONVERT)

$(LBL2NC): $(LBL2NC_OBJECTS)
	$(FC) $(FCFLAGS) $(LBL2NC_OBJECTS) -o $(LBL2NC) $(NETCDF_LIB) $(OMPFLAG)

$(CONTINUUM): $(CONTINUUM_OBJECTS)
	$(FC) $(FCFLAGS) $(CONTINUUM_OBJECTS) -o $(CONTINUUM) $(NETCDF_LIB)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

interpolation.o: parkind1.o
ckdmip_tool.o: parkind1.o easy_netcdf.o interpolation.o
ckdmip_convert.o: parkind1.o easy_netcdf.o
ckdmip_sw.o: shortwave_fluxes.o easy_netcdf.o shortwave_config.o ckdmip_version.h
shortwave_fluxes.o: parkind1.o
change_continuum.o: parkind1.o easy_netcdf.o caviar_continuum.o interpolation.o continuum_lut.o
continuum_lut.o: parkind1.o interpolation.o

ckdmip_lw.o: longwave_config.o parkind1.o easy_netcdf.o planck_function.o \
	longwave_fluxes.o ckdmip_version.h
lblrtm2nc.o: parkind1.o easy_netcdf.o
longwave_fluxes.o: parkind1.o
planck_function.o: parkind1.o
longwave_config.o: parkind1.o
easy_netcdf.o: parkind1.o radiation_io.o
radiation_io.o: yomlun_ifsaux.o abor1.o
abor1.o: parkind1.o yomlun_ifsaux.o

clean:
	rm -f $(OBJECTS) *.mod $(PROGRAMS)
