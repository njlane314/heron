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
IO_SRC = framework/io/src/ArtFileProvenanceIO.cpp \
         framework/io/src/EventListIO.cpp \
         framework/io/src/NormalisationService.cpp \
         framework/io/src/RunDatabaseService.cpp \
         framework/io/src/SnapshotService.cpp \
         framework/io/src/SampleIO.cpp \
         framework/io/src/SubRunInventoryService.cpp
OBJ_DIR = build/obj
IO_OBJ = $(IO_SRC:%.cpp=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = build/lib/libheronAna.so
ANA_SRC = framework/ana/src/AnalysisModel.cpp \
          framework/ana/src/AnalysisConfigService.cpp \
          framework/ana/src/ColumnDerivationService.cpp \
          framework/ana/src/EventSampleFilterService.cpp \
          framework/ana/src/RDataFrameService.cpp \
          framework/ana/src/SelectionService.cpp
ANA_OBJ = $(ANA_SRC:%.cpp=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = build/lib/libheronPlot.so
PLOT_SRC = framework/plot/src/Plotter.cpp \
           framework/plot/src/StackedHist.cpp \
           framework/plot/src/UnstackedHist.cpp \
           framework/plot/src/PlottingHelper.cpp \
           framework/plot/src/EfficiencyPlot.cpp \
           framework/plot/src/TemplateBinningOptimiser1D.cpp
PLOT_OBJ = $(PLOT_SRC:%.cpp=$(OBJ_DIR)/%.o)

EVD_LIB_NAME = build/lib/libheronEvd.so
EVD_SRC = framework/evd/src/EventDisplay.cpp
EVD_OBJ = $(EVD_SRC:%.cpp=$(OBJ_DIR)/%.o)

HERON_NAME = build/bin/heron
APPS_SRC = framework/apps/src/heron.cpp \
           framework/apps/src/ArtWorkflow.cpp \
           framework/apps/src/SampleWorkflow.cpp \
           framework/apps/src/EventWorkflow.cpp
APPS_OBJ = $(APPS_SRC:%.cpp=$(OBJ_DIR)/%.o)

INCLUDES = -I./framework/core/include -I./framework/io/include -I./framework/ana/include -I./framework/plot/include -I./framework/evd/include -I./framework/apps/include

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

$(HERON_NAME): $(APPS_OBJ) $(IO_LIB_NAME) $(ANA_LIB_NAME) $(PLOT_LIB_NAME)
	mkdir -p $(dir $(HERON_NAME))
	$(CXX) $(CXXFLAGS) $(APPS_OBJ) -Lbuild/lib -lheronIO \
		-lheronAna -lheronPlot $(LDFLAGS) -o $(HERON_NAME)

$(OBJ_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -c $< -o $@

clean:
	rm -rf build/lib build/bin build/obj
