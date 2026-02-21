CXX ?= $(shell $(ROOT_CONFIG) --cxx)
ROOT_CONFIG ?= root-config

FRAMEWORK_DIR ?= framework
BUILD_DIR ?= build
LIB_DIR := $(BUILD_DIR)/lib
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj

ifeq ($(shell command -v $(ROOT_CONFIG) 2>/dev/null),)
$(error ROOT not found. Please set up ROOT so that 'root-config' is on your PATH.)
endif

NLOHMANN_JSON_INC ?=
NLOHMANN_JSON_CFLAGS ?=
ifneq ($(strip $(NLOHMANN_JSON_INC)),)
NLOHMANN_JSON_CFLAGS := -isystem $(NLOHMANN_JSON_INC)
endif

CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra $(shell $(ROOT_CONFIG) --cflags) $(NLOHMANN_JSON_CFLAGS)
LDFLAGS ?= $(shell $(ROOT_CONFIG) --libs) -lsqlite3

MODULES := io ana plot evd
MODULES_DIR := $(FRAMEWORK_DIR)/modules

IO_LIB_NAME = $(LIB_DIR)/libheronIO.so
IO_SRC = $(MODULES_DIR)/io/src/ArtFileProvenanceIO.cc \
         $(MODULES_DIR)/io/src/EventListIO.cc \
         $(MODULES_DIR)/io/src/NormalisationService.cc \
         $(MODULES_DIR)/io/src/RunDatabaseService.cc \
         $(MODULES_DIR)/io/src/SnapshotService.cc \
         $(MODULES_DIR)/io/src/SampleIO.cc \
         $(MODULES_DIR)/io/src/SubRunInventoryService.cc
IO_OBJ = $(IO_SRC:%.cc=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = $(LIB_DIR)/libheronAna.so
ANA_SRC = $(MODULES_DIR)/ana/src/AnalysisConfigService.cc \
          $(MODULES_DIR)/ana/src/ColumnDerivationService.cc \
          $(MODULES_DIR)/ana/src/EventSampleFilterService.cc \
          $(MODULES_DIR)/ana/src/RDataFrameService.cc \
          $(MODULES_DIR)/ana/src/SelectionService.cc
ANA_OBJ = $(ANA_SRC:%.cc=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = $(LIB_DIR)/libheronPlot.so
PLOT_SRC = $(MODULES_DIR)/plot/src/Plotter.cc \
           $(MODULES_DIR)/plot/src/StackedHist.cc \
           $(MODULES_DIR)/plot/src/UnstackedHist.cc \
           $(MODULES_DIR)/plot/src/PlottingHelper.cc \
           $(MODULES_DIR)/plot/src/EfficiencyPlot.cc
PLOT_OBJ = $(PLOT_SRC:%.cc=$(OBJ_DIR)/%.o)

EVD_LIB_NAME = $(LIB_DIR)/libheronEvd.so
EVD_SRC = $(MODULES_DIR)/evd/src/EventDisplay.cc
EVD_OBJ = $(EVD_SRC:%.cc=$(OBJ_DIR)/%.o)

HERON_NAME = $(BIN_DIR)/heron
CORE_SRC = $(FRAMEWORK_DIR)/core/src/heron.cc \
           $(FRAMEWORK_DIR)/core/src/ArtWorkflow.cc \
           $(FRAMEWORK_DIR)/core/src/SampleWorkflow.cc \
           $(FRAMEWORK_DIR)/core/src/EventWorkflow.cc
CORE_OBJ = $(CORE_SRC:%.cc=$(OBJ_DIR)/%.o)

INCLUDES = $(addprefix -I./$(MODULES_DIR)/,$(addsuffix /include,$(MODULES))) -I./$(FRAMEWORK_DIR)/core/include

define build_shared_library
$(1): $$($(2)_OBJ)
	mkdir -p $$(dir $$@)
	$$(CXX) -shared $$(CXXFLAGS) $$^ $$(LDFLAGS) -o $$@
endef

all: $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME) $(EVD_LIB_NAME) $(HERON_NAME)

$(eval $(call build_shared_library,$(IO_LIB_NAME),IO))
$(eval $(call build_shared_library,$(ANA_LIB_NAME),ANA))
$(eval $(call build_shared_library,$(PLOT_LIB_NAME),PLOT))
$(eval $(call build_shared_library,$(EVD_LIB_NAME),EVD))

$(HERON_NAME): $(CORE_OBJ) $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME)
	mkdir -p $(dir $(HERON_NAME))
	$(CXX) $(CXXFLAGS) $(CORE_OBJ) -L$(LIB_DIR) -lheronIO \
		-lheronAna -lheronPlot $(LDFLAGS) -o $(HERON_NAME)

$(OBJ_DIR)/%.o: %.cc
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf $(LIB_DIR) $(BIN_DIR) $(OBJ_DIR)
