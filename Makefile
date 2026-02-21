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
IO_SRC = io/src/ArtFileProvenanceIO.cpp \
         io/src/EventListIO.cpp \
         io/src/NormalisationService.cpp \
         io/src/RunDatabaseService.cpp \
         io/src/SnapshotService.cpp \
         io/src/SampleIO.cpp \
         io/src/SubRunInventoryService.cpp
OBJ_DIR = build/obj
IO_OBJ = $(IO_SRC:%.cpp=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = build/lib/libheronAna.so
ANA_SRC = ana/src/AnalysisConfigService.cpp \
          ana/src/ColumnDerivationService.cpp \
          ana/src/EventSampleFilterService.cpp \
          ana/src/RDataFrameService.cpp \
          ana/src/SelectionService.cpp
ANA_OBJ = $(ANA_SRC:%.cpp=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = build/lib/libheronPlot.so
PLOT_SRC = plot/src/Plotter.cpp \
           plot/src/StackedHist.cpp \
           plot/src/UnstackedHist.cpp \
           plot/src/PlottingHelper.cpp \
           plot/src/EfficiencyPlot.cpp \
           plot/src/TemplateBinningOptimiser1D.cpp
PLOT_OBJ = $(PLOT_SRC:%.cpp=$(OBJ_DIR)/%.o)

EVD_LIB_NAME = build/lib/libheronEvd.so
EVD_SRC = evd/src/EventDisplay.cpp
EVD_OBJ = $(EVD_SRC:%.cpp=$(OBJ_DIR)/%.o)

HERON_NAME = build/bin/heron
APPS_SRC = apps/src/heron.cpp \
           apps/src/ArtWorkflow.cpp \
           apps/src/SampleWorkflow.cpp \
           apps/src/EventWorkflow.cpp
APPS_OBJ = $(APPS_SRC:%.cpp=$(OBJ_DIR)/%.o)

INCLUDES = -I./io/include -I./ana/include -I./plot/include -I./evd/include -I./apps/include

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
