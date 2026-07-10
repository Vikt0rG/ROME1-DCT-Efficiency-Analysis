INCLUDE_DIR := include
SRC_DIR     := src
UTILS_DIR   := utils
DATA_DIR    := data
BIN_DIR     := bin
OBJ_DIR     := obj

EXECUTABLE_NAME := analysis
TEST_DATA_DIR   := $(DATA_DIR)/example/input/

# Map out all sources matching your exact new nested layout
SRCS := $(SRC_DIR)/main.cpp \
        $(SRC_DIR)/core/hit.cpp \
        $(SRC_DIR)/core/cluster.cpp \
        $(SRC_DIR)/core/track.cpp \
        $(SRC_DIR)/core/event.cpp \
        $(SRC_DIR)/analysis/dataProcesser.cpp \
        $(SRC_DIR)/analysis/dataAnalyzer.cpp \
        $(SRC_DIR)/plotting/dataPlotter.cpp \
        $(SRC_DIR)/plotting/plotStyler.cpp \
        $(SRC_DIR)/plotting/plotBatchExporter.cpp

# Convert source paths to flat object targets (e.g., src/core/hit.cpp -> obj/hit.o)
OBJS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRCS)))

# Compiler configuration via ROOT settings
CC := $(shell root-config --cxx)

CONDA_PREFIX_CLEAN = $(strip $(CONDA_PREFIX))
EXTRA_INCLUDES :=
EXTRA_LIBS :=

ifneq ($(CONDA_PREFIX_CLEAN),)
    EXTRA_INCLUDES := -I$(CONDA_PREFIX_CLEAN)/include
    EXTRA_LIBS     := -L$(CONDA_PREFIX_CLEAN)/lib
    LDFLAGS        += $(EXTRA_LIBS)
endif

CFLAGS  := -g -Wall -std=c++17 -O2 \
           -I$(INCLUDE_DIR) \
           -I$(INCLUDE_DIR)/core \
           -I$(INCLUDE_DIR)/analysis \
           -I$(INCLUDE_DIR)/plotting \
           -I$(UTILS_DIR) \
           $(EXTRA_INCLUDES) \
           $(shell root-config --cflags)

LDFLAGS := $(sort $(shell root-config --libs) -lyaml-cpp)

# Tell Make where to search for raw source files dynamically across your subsystems
vpath %.cpp $(SRC_DIR) $(SRC_DIR)/core $(SRC_DIR)/analysis $(SRC_DIR)/plotting

# Default targets
all: build

# Link Step
build: $(BIN_DIR)/$(EXECUTABLE_NAME)

$(BIN_DIR)/$(EXECUTABLE_NAME): $(OBJS)
	@echo ""
	@echo "Linking final executable..."
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)
	@echo ""
	@echo "Build complete. Executable is at $@"
	@echo ""

# Compilation Pattern Rule for flattening objects into /obj
$(OBJ_DIR)/%.o: %.cpp
	@echo "Compiling $<..."
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo ""
	@echo "Cleaning build artifacts..."
	rm -rf $(BIN_DIR) $(OBJ_DIR) output.root
	@echo "Clean complete."
	@echo ""

rebuild: clean build

help:
	@echo "ROME1-DCT-Efficiency-Analysis Build System"
	@echo "===================================="
	@echo "Available targets:"
	@echo "  build          - Compile the executable (default)"
	@echo "  clean          - Remove build artifacts"
	@echo "  help           - Show this help message"

.PHONY: all build clean help