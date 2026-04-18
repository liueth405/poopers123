################################################################################
######################### User configurable parameters #########################
# filename extensions
CEXTS:=c
ASMEXTS:=s S
CXXEXTS:=cpp c++ cc

# probably shouldn't modify these, but you may need them below
ROOT=.
FWDIR:=$(ROOT)/firmware
BINDIR=$(ROOT)/bin
SRCDIR=$(ROOT)/src
COLDSRCDIR=$(ROOT)/cold_src
INCDIR=$(ROOT)/include

WARNFLAGS+=
EXTRA_CFLAGS=-O3
EXTRA_CXXFLAGS=-O3

# Add third-party include directories to the search path
# This allows using <blasfeo_target.h> instead of "blasfeo/blasfeo_target.h"
# and is required for third-party headers that use angle brackets to find each other.
THIRD_PARTY_INCLUDES = -I$(INCDIR)/blasfeo -I$(INCDIR)/hpipm/include -I$(INCDIR)/acados -I$(INCDIR)/acados_c -I$(INCDIR)/mpcc
EXTRA_CFLAGS += $(THIRD_PARTY_INCLUDES)
EXTRA_CXXFLAGS += $(THIRD_PARTY_INCLUDES)

# Set to 1 to enable hot/cold linking
USE_PACKAGE:=1

# Add libraries you do not wish to include in the cold image here
# EXCLUDE_COLD_LIBRARIES:= $(FWDIR)/your_library.a
EXCLUDE_COLD_LIBRARIES:=
# Add our custom cold math library to the links
LIBAR = $(BINDIR)/$(LIBNAME).a
LIBRARIES += $(LIBAR)
COLD_LIBRARIES += $(LIBAR) $(BINDIR)/cold_math.cpp.o

# Set this to 1 to add additional rules to compile your project as a PROS library template
IS_LIBRARY:=0
# TODO: CHANGE THIS!
# Be sure that your header files are in the include directory inside of a folder with the
# same name as what you set LIBNAME to below.
LIBNAME:=libbest
VERSION:=1.0.0
# EXCLUDE_SRC_FROM_LIB= $(SRCDIR)/unpublishedfile.c
# this line excludes opcontrol.c and similar files
EXCLUDE_SRC_FROM_LIB+=$(foreach file, $(SRCDIR)/main,$(foreach cext,$(CEXTS),$(file).$(cext)) $(foreach cxxext,$(CXXEXTS),$(file).$(cxxext)))

# files that get distributed to every user (beyond your source archive) - add
# whatever files you want here. This line is configured to add all header files
# that are in the directory include/LIBNAME
TEMPLATE_FILES=$(INCDIR)/$(LIBNAME)/*.h $(INCDIR)/$(LIBNAME)/*.hpp

# Rule for Hot Image Objects (only from src/)
ELF_DEPS = $(sort $(foreach cext,$(CEXTS),$(call rwildcard,$(SRCDIR),*.$(cext))) \
            $(foreach cxxext,$(CXXEXTS),$(call rwildcard,$(SRCDIR),*.$(cxxext))))
ELF_DEPS := $(addprefix $(BINDIR)/,$(patsubst $(SRCDIR)/%,%.o,$(ELF_DEPS)))

# Rule for Cold Library Objects (everything in src/ except main and MPCC related, plus cold_src)
# First get standard src/ objects using standard PROS wildcards:
# We EXCLUDE the files we are actively tuning to keep them in the HOT package for fast uploads.
HOT_FILES = $(BINDIR)/main.cpp.o \
	$(BINDIR)/motions/motion.cpp.o \
	$(BINDIR)/motions/mpcc_controller.cpp.o \
	$(BINDIR)/motions/ilqr_controller.cpp.o \
	$(BINDIR)/motions/mpcc_sd_log.cpp.o \
	$(BINDIR)/motions/pathplanner.cpp.o

LIB_OBJS = $(filter-out $(HOT_FILES), $(ELF_DEPS))
# Add our custom cold code explicitly:
LIB_OBJS += $(BINDIR)/cold_math.cpp.o

.DEFAULT_GOAL=quick

########## Nothing below this line should be edited by typical users ###########
-include ./common.mk

# Custom rules (must be at the bottom to override common.mk)
$(BINDIR)/%.o: $(COLDSRCDIR)/%
	$(VV)mkdir -p $(dir $@)
	$(call test_output_2,Compiled Cold $$< ,$(CXX) -c $(INCLUDE) $(CXXFLAGS) -o $@ $<,$(OK_STRING))

# Override the library rule to use our specific LIB_OBJS
$(LIBAR): $(LIB_OBJS)
	-$Dmkdir $(BINDIR)
	-$Drm -f $@
	$(call test_output_2,Creating $@ ,$(AR) rcs $@ $^, $(DONE_STRING))
