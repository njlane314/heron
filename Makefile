CXX ?= $(shell $(ROOT_CONFIG) --cxx)
ROOT_CONFIG ?= root-config

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

IO_LIB_NAME = build/lib/libheronIO.so
IO_SRC = framework/io/src/ArtFileProvenanceIO.cc \
         framework/io/src/EventListIO.cc \
         framework/io/src/NormalisationService.cc \
         framework/io/src/RunDatabaseService.cc \
         framework/io/src/SnapshotService.cc \
         framework/io/src/SampleIO.cc \
         framework/io/src/SubRunInventoryService.cc
OBJ_DIR = build/obj
IO_OBJ = $(IO_SRC:%.cc=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = build/lib/libheronAna.so
ANA_SRC = framework/ana/src/AnalysisConfigService.cc \
          framework/ana/src/ColumnDerivationService.cc \
          framework/ana/src/EventSampleFilterService.cc \
          framework/ana/src/RDataFrameService.cc \
          framework/ana/src/SelectionService.cc
ANA_OBJ = $(ANA_SRC:%.cc=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = build/lib/libheronPlot.so
PLOT_SRC = framework/plot/src/Plotter.cc \
           framework/plot/src/StackedHist.cc \
           framework/plot/src/UnstackedHist.cc \
           framework/plot/src/PlottingHelper.cc \
           framework/plot/src/EfficiencyPlot.cc
PLOT_OBJ = $(PLOT_SRC:%.cc=$(OBJ_DIR)/%.o)

EVD_LIB_NAME = build/lib/libheronEvd.so
EVD_SRC = framework/evd/src/EventDisplay.cc
EVD_OBJ = $(EVD_SRC:%.cc=$(OBJ_DIR)/%.o)

HERON_NAME = build/bin/heron
CORE_SRC = framework/core/src/heron.cc \
           framework/core/src/ArtWorkflow.cc \
           framework/core/src/SampleWorkflow.cc \
           framework/core/src/EventWorkflow.cc
CORE_OBJ = $(CORE_SRC:%.cc=$(OBJ_DIR)/%.o)

INCLUDES = -I./framework/io/include -I./framework/ana/include -I./framework/plot/include -I./framework/evd/include -I./framework/core/include

all: $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME) $(EVD_LIB_NAME) $(HERON_NAME)

$(IO_LIB_NAME): $(IO_OBJ)
	mkdir -p $(dir $(IO_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(IO_OBJ) $(LDFLAGS) -o $(IO_LIB_NAME)

$(ANA_LIB_NAME): $(ANA_OBJ)
	mkdir -p $(dir $(ANA_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(ANA_OBJ) $(LDFLAGS) -o $(ANA_LIB_NAME)

$(PLOT_LIB_NAME): $(PLOT_OBJ)
	mkdir -p $(dir $(PLOT_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(PLOT_OBJ) $(LDFLAGS) -o $(PLOT_LIB_NAME)

$(EVD_LIB_NAME): $(EVD_OBJ)
	mkdir -p $(dir $(EVD_LIB_NAME))
	$(CXX) -shared $(CXXFLAGS) $(EVD_OBJ) $(LDFLAGS) -o $(EVD_LIB_NAME)

$(HERON_NAME): $(CORE_OBJ) $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME)
	mkdir -p $(dir $(HERON_NAME))
	$(CXX) $(CXXFLAGS) $(CORE_OBJ) -Lbuild/lib -lheronIO \
		-lheronAna -lheronPlot $(LDFLAGS) -o $(HERON_NAME)

$(OBJ_DIR)/%.o: %.cc
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build/lib build/bin build/obj
