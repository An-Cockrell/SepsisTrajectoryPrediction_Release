src := agentConstructors.cpp cellFunctions.cpp IIRABM_ctypes.cpp simulationFunctions.cpp run.cpp

CXX = mpicxx
CXXLD = mpicxx

SWIG = /usr/bin/swig
SED = sed

CXX_RELEASE_FLAGS = -Wall -O2 -g0 -std=c++11 -MMD -MP
CXX_DEBUG_FLAGS = -Wall -O0 -g3 -std=c++11 -MMD -MP

CXX_FLAGS = $(CXX_RELEASE_FLAGS)

# these are not necessary unless you are creating a tcl wrapper
# to use with swift/T
TCL_INCLUDE =  /usr/include/tcl8.6
TCL_LIB     = -L/usr/lib/x86_64-linux-gnu -ltcl8.6
TCLSH       = /usr/bin/tclsh8.6
TCLSH_LOCAL = /usr/bin/tclsh8.6

INCLUDES := -I ../
INCLUDES += -I $(TCL_INCLUDE)

LIB_EXT =

OS := $(shell uname)

ifeq ($(OS),Darwin)
 LIB_LDFLAGS += -dynamiclib
 LIB_EXT = dylib
else
 # linux
 LIB_CXX_FLAGS += -fPIC
 LIB_LDFLAGS += -shared 
 LIB_EXT = so
endif

SRC_DIR= ../
BUILD_DIR = ./build

# objects used by both executable and tests
OBJECTS := $(subst .cpp,.o, $(addprefix $(BUILD_DIR)/, $(src)))

VPATH = ../

NAME = iirabm
TCL_LIB_NAME=libtcl$(NAME).$(LIB_EXT)
LIB_NAME=$(NAME).$(LIB_EXT)

SED := sed
MV := mv -f

SWIFT_MODEL_DIR=../../swift_proj/model

-include $(OBJECTS:.o=.d)

.PHONY: all clean

all: iirabm_tcl_lib

iirabm_lib: $(OBJECTS)
	$(CXXLD)  $(filter-out %.d, $^) $(LIB_LDFLAGS)  -o ./$(LIB_NAME)

$(BUILD_DIR)/%.o : %.cpp
	@-mkdir -p $(dir $@)
	$(CXX) $(CXX_FLAGS) $(LIB_CXX_FLAGS) $(INCLUDES)  -c $< -o $@

iirabm_tcl_lib: $(OBJECTS) $(BUILD_DIR)/iirabm_wrapper.o
	$(CXXLD)  $(filter-out %.d, $^) $(LIB_LDFLAGS) $(LIBS) $(TCL_LIB) -o $(SWIFT_MODEL_DIR)/$(TCL_LIB_NAME)
	LEAF_PKG=iirabm_model LEAF_VERSION=0.0 \
	LEAF_SO=$(TCL_LIB_NAME) LEAF_TCL=iirabm_model.tcl \
	$(TCLSH_LOCAL) $(SWIFT_MODEL_DIR)/make-package.tcl > $(SWIFT_MODEL_DIR)/pkgIndex.tcl
	
../iirabm_wrapper.cpp :  $(SWIFT_MODEL_DIR)/iirabm_model.i ./run.h
	$(SWIG) -c++ -o $(@) $(<)
	$(SED) -i.bu "s/Iirabm_model_Init/Tcliirabm_Init/" $(@)
		

clean:
	rm -fv $(LIB_NAME)
	rm -rf $(BUILD_DIR)/*
	rm -rf ../swift/lib/$(TCL_LIB_NAME)
