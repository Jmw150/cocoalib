# Makefile configuration for CoCoALib.
# Created automatically by the configure script.
# Created on  2022-03-23  at time  15:39:35
# Command was: 
# ./configure 

PLATFORM=Linux 5.13.0-37-generic x86_64

######################################################
# Definitions common to all Makefiles in CoCoALib.

# Version number of CoCoALib we shall build.
COCOALIB_VERSION=0.99800

INSTALL_CMD=install
COCOALIB_INSTALL_DIR=/usr/local
COCOA5_INSTALL_DIR=/usr/local/bin

EXTLIBS=$(COCOA_ROOT)/configuration/ExternalLibs


######################################################
# Compilation settings.

CXX=g++
CXXFLAGS=  -std=c++14  -Wall  -pedantic  -fPIC  -O2  
# for debugging append flags -Og  -g   (-pg profiling)

ARFLAGS=-rcS

ISYSTEM=isystem

######################################################
# External libraries for CoCoALib:

# We use the following GMP installation:
GMP_VERSION=6.2.0
GMP_LDLIB=-lgmp
GMPXX_LDLIB=-lgmpxx

# FROBBY settings:
HAVE_FROBBY=no

# CDD settings:
HAVE_CDD=no

# GFAN settings:
HAVE_GFAN=no

# GSL settings:
HAVE_GSL=no

# MathSAT settings:
HAVE_MATHSAT =no

# Normaliz settings:
HAVE_NORMALIZ=no

LDLIBS= $(COCOA_LIB) -L$(EXTLIBS)/lib  $(FROBBY_LDLIBS)  $(GFAN_LDLIBS)  $(CDD_LDLIBS)  $(GSL_LDLIBS)  $(MATHSAT_LDLIBS)  $(NORMALIZ_LDLIBS)  $(GMPXX_LDLIB)  $(GMP_LDLIB)


########################################
# Flags/libraries/settings for CoCoA-5

# BOOST library:
HAVE_BOOST=yes
BOOST_LDLIBS=-lboost_filesystem-symlink -lboost_thread-symlink -lboost_system-symlink 

# READLINE settings:
HAVE_READLINE=yes
READLINE_LDLIBS=-lreadline

COCOA5_CXX_DEFINES=-DCoCoA_WITH_READLINE
COCOA5_LDLIBS=$(BOOST_LDLIBS)  $(READLINE_LDLIBS)

# Flag saying whether to build the Qt GUI
BUILD_QT_GUI=no



##################################################################
# This is the second fixed part of the common Makefile definitions.
# The configure script will copy it to autoconf.mk.
# NOTE: COCOA_ROOT is defined as a relative path in each individual Makefile.

COCOA_HDR=$(COCOA_ROOT)/include/CoCoA/library.H

INCLUDE_PATHS=-I$(COCOA_ROOT)/include  -$(ISYSTEM) $(EXTLIBS)/include
COMPILE=$(CXX)  $(CXXFLAGS)  $(INCLUDE_PATHS)
COCOA_LIB=$(COCOA_ROOT)/lib/libcocoa.a



# Rule for compiling C++ code in *.C files into *.o object files
%.o: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE) -c -o $@ $<

# Rule for compiling and linking C++ code in *.C files
%: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE) -o $@ $< $(LDLIBS)
	@AppleDir="$@.dSYM" ; \
	echo " " $(CXXFLAGS) " " | fgrep " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


# Rule for compiling C++ code in *.cpp files into *.o object files
%.o: %.cpp
	@echo "Compiling `basename $@`"
	$(COMPILE) -c -o $@ $<

# Rule for compiling and linking C++ code in *.cpp files
%: %.cpp
	@echo "Compiling `basename $@`"
	$(COMPILE) -o $@ $< $(LDLIBS)
	@AppleDir="$@.dSYM" ; \
	echo " " $(CXXFLAGS) " " | fgrep " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


# The following are derived from the conventions for Makefiles for code
# which should become part of the GNU project.  It seems reasonable to
# adopt them here as well.  I found the recommendations in the online
# help for make (e.g. run the command "info make" on a GNU/Linux system).
SHELL=/bin/bash
.SUFFIXES:
.SUFFIXES: .C .o
