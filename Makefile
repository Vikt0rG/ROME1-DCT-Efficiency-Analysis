INCLUDE_DIR := include
SRC_DIR := src
UTILS_DIR := utils
DATA_DIR := data
BIN_DIR := bin
OBJ_DIR := obj

EXECUTABLE_NAME := analysis
TEST_DATA_DIR := $(DATA_DIR)/example/input/

# Collect all source files from src/
SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.cxx)

# Collect all header files from include/
HEADERS := $(wildcard $(INCLUDE_DIR)/*.hpp)

# Compiler: use the same C++ compiler that ROOT was built with
CC = $(shell root-config --cxx)

CONDA_PREFIX_CLEAN = $(strip $(CONDA_PREFIX))
EXTRA_INCLUDES :=
EXTRA_LIBS :=

ifneq ($(CONDA_PREFIX_CLEAN),)
    EXTRA_INCLUDES := -I$(CONDA_PREFIX_CLEAN)/include
    EXTRA_LIBS := -L$(CONDA_PREFIX_CLEAN)/lib -Wl,-rpath,/$(CONDA_PREFIX_CLEAN)/lib
endif

# Compiler flags:
# -g 						: debugging info
# -Wall 					: compiler warnings
# -O2 						: optimization level 2
# -std=c++17 				: use C++ standard as the one from root-config
# `root-config --cflags` 	: necessary ROOT include flags
CFLAGS := -g -Wall -std=c++17 -O2 -I$(INCLUDE_DIR) -I$(UTILS_DIR) $(EXTRA_INCLUDES) $(shell root-config --cflags)

LDFLAGS := $(shell root-config --libs) -lHistPainter -lyaml-cpp

# Default target
all: build

build: $(SRCS) $(HEADERS)
	@echo "Compiling DCT analysis source files..."
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(SRCS) -o $(BIN_DIR)/$(EXECUTABLE_NAME) $(LDFLAGS)
	@echo "Build complete. Executable is at $(BIN_DIR)/$(EXECUTABLE_NAME)"

test: build
	@echo "Running DCT analysis with test data..."
	@echo "Input directory: $(TEST_DATA_DIR)"
	@echo "Parameters: dt_max=-100, dt_min=-180"
	./$(BIN_DIR)/$(EXECUTABLE_NAME) $(TEST_DATA_DIR) -100 -180

clean:
	@echo "Cleaning build artifacts..."
	@if [ -d $(BIN_DIR) ]; then rm -rf $(BIN_DIR); fi
	@if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@rm -f output.root
	@echo "Clean complete."

help:
	@echo "DCT Efficiency Analysis Build System"
	@echo "===================================="
	@echo "Available targets:"
	@echo "  build          - Compile the executable (default)"
	@echo "  test           - Build and run with test data"
	@echo "  clean          - Remove build artifacts"
	@echo "  help           - Show this help message"

.PHONY: all build test clean help