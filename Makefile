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

all: build process_raw_hits

build: $(SRCS)
	@echo "Compiling analysis source files..."
	@mkdir -p bin
	$(CC) $(CFLAGS) $(filter-out $(SRC_DIR)/main.cpp $(SRC_DIR)/process_raw_hits_utils.cpp, $(SRCS)) -o bin/$(EXECUTABLE_NAME) $(shell root-config --libs)
	@echo "Build complete. Executable is at bin/$(EXECUTABLE_NAME)"

process_raw_hits: $(SRC_DIR)/main.cpp $(SRC_DIR)/process_raw_hits_utils.cpp
	@echo "Compiling process_raw_hits standalone executable..."
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $^ -o $(BIN_DIR)/$(PROCESS_RAW_HITS) $(shell root-config --libs)
	@echo "Build complete. Executable is at $(BIN_DIR)/$(PROCESS_RAW_HITS)"

clean:
	@echo "Cleaning build artifacts..."
	@if [ -d bin ]; then rm -rf bin; fi
	@if [ -d obj ]; then rm -rf obj; fi
	@echo "Clean complete."