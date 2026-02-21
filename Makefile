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
IO_SRC = mod/io/src/ArtFileProvenanceIO.cpp \
         mod/io/src/EventListIO.cpp \
         mod/io/src/NormalisationService.cpp \
         mod/io/src/RunDatabaseService.cpp \
         mod/io/src/SnapshotService.cpp \
         mod/io/src/SampleIO.cpp \
         mod/io/src/SubRunInventoryService.cpp
OBJ_DIR = build/obj
IO_OBJ = $(IO_SRC:%.cpp=$(OBJ_DIR)/%.o)

ANA_LIB_NAME = build/lib/libheronAna.so
ANA_SRC = mod/ana/src/AnalysisConfigService.cpp \
          mod/ana/src/ColumnDerivationService.cpp \
          mod/ana/src/EventSampleFilterService.cpp \
          mod/ana/src/RDataFrameService.cpp \
          mod/ana/src/SelectionService.cpp
ANA_OBJ = $(ANA_SRC:%.cpp=$(OBJ_DIR)/%.o)

PLOT_LIB_NAME = build/lib/libheronPlot.so
PLOT_SRC = mod/plot/src/Plotter.cpp \
           mod/plot/src/StackedHist.cpp \
           mod/plot/src/UnstackedHist.cpp \
           mod/plot/src/PlottingHelper.cpp \
           mod/plot/src/EfficiencyPlot.cpp \
           mod/plot/src/TemplateBinningOptimiser1D.cpp
PLOT_OBJ = $(PLOT_SRC:%.cpp=$(OBJ_DIR)/%.o)

EVD_LIB_NAME = build/lib/libheronEvd.so
EVD_SRC = mod/evd/src/EventDisplay.cpp
EVD_OBJ = $(EVD_SRC:%.cpp=$(OBJ_DIR)/%.o)

HERON_NAME = build/bin/heron
APPS_SRC = mod/apps/src/heron.cpp \
           mod/apps/src/ArtWorkflow.cpp \
           mod/apps/src/SampleWorkflow.cpp \
           mod/apps/src/EventWorkflow.cpp
APPS_OBJ = $(APPS_SRC:%.cpp=$(OBJ_DIR)/%.o)

INCLUDES = -I./mod/io/include -I./mod/ana/include -I./mod/plot/include -I./mod/evd/include -I./mod/apps/include

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
