INCLUDE_DIR := include/
SRC_DIR := src/
BIN_DIR := bin
OBJ_DIR := obj

EXECUTABLE_NAME := analysis
PROCESS_RAW_HITS := process_raw_hits

# Collect all source files from src/
SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.cxx)

# Compiler: use the same C++ compiler that ROOT was built with
CC = $(shell root-config --cxx)

# Compiler flags:
# -g 						: debugging info
# -Wall 					: compiler warnings
# -O2 						: optimization level 2
# -std=c++17 				: use C++ standard as the one from root-config
# `root-config --cflags` 	: necessary ROOT include flags
CFLAGS := -g -Wall -std=c++17 -O2 -I$(INCLUDE_DIR) -I$(UTIL_DIR) $(shell root-config --cflags)

LDFLAGS :=
all: build process_raw_hits

build: $(SRCS)
	@echo "Compiling analysis source files..."
	@mkdir -p bin
	$(CC) $(CFLAGS) $(filter-out $(SRC_DIR)/main.cpp $(SRC_DIR)/process_raw_hits_utils.cpp, $(SRCS)) -o bin/$(EXECUTABLE_NAME) $(shell root-config --libs) $(LDFLAGS)
	@echo "Build complete. Executable is at bin/$(EXECUTABLE_NAME)"

process_raw_hits: $(SRC_DIR)/main.cpp $(SRC_DIR)/process_raw_hits_utils.cpp
	@echo "Compiling process_raw_hits standalone executable..."
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $^ -o $(BIN_DIR)/$(PROCESS_RAW_HITS) $(shell root-config --libs) $(LDFLAGS)
	@echo "Build complete. Executable is at $(BIN_DIR)/$(PROCESS_RAW_HITS)"

clean:
	@echo "Cleaning build artifacts..."
	@if [ -d bin ]; then rm -rf bin; fi
	@if [ -d obj ]; then rm -rf obj; fi
	@echo "Clean complete."

# Draft directory compilation
DRAFT_DIR := draft
DRAFT_SRC_DIR := $(DRAFT_DIR)/src
DRAFT_INCLUDE_DIR := $(DRAFT_DIR)/include
DRAFT_UTILS_DIR := $(DRAFT_DIR)/utils
DRAFT_DATA_DIR := $(DRAFT_DIR)/data
DRAFT_EXECUTABLE_NAME := process_raw_hits_draft

# Collect all source files from draft/src/
DRAFT_SRCS := $(wildcard $(DRAFT_SRC_DIR)/*.cpp)

draft_build: $(DRAFT_SRCS)
	@echo "Compiling draft source files..."
	@mkdir -p bin
	$(CC) $(CFLAGS) -I$(DRAFT_INCLUDE_DIR) -I$(DRAFT_DIR) $(DRAFT_SRCS) -o bin/$(DRAFT_EXECUTABLE_NAME) $(shell root-config --libs) $(LDFLAGS)
	@echo "Build complete. Executable is at bin/$(DRAFT_EXECUTABLE_NAME)"

draft_test: draft_build
	@echo "Running draft test with parameters: input=$(DRAFT_DATA_DIR) dt_max=-100 dt_min=-180"
	./bin/$(DRAFT_EXECUTABLE_NAME) $(DRAFT_DATA_DIR) -100 -180

draft_clean:
	@echo "Cleaning draft build artifacts..."
	@rm -f bin/$(DRAFT_EXECUTABLE_NAME)
	@echo "Draft clean complete."

.PHONY: all build process_raw_hits clean draft_build draft_test draft_clean