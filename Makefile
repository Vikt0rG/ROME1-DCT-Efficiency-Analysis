INCLUDE_DIR := include/
SRC_DIR := src/

EXECUTABLE_NAME := analysis

# Collect all source files from src/
SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.cxx)

# Compiler: use the same C++ compiler that ROOT was built with
CC = $(shell root-config --cxx)

# Compiler flags:
# -g 						: enables debugging information
# -Wall 					: turns on most, but not all, compiler warnings
# -O0 						: no optimization (easier debugging)
# -std=c++17 				: use C++ standard as the one from root-config
# `root-config --cflags` 	: adds necessary ROOT include flags
CFLAGS := -g -Wall -std=c++17 -O0 -I$(INCLUDE_DIR) -I$(UTIL_DIR) $(shell root-config --cflags)

all: build

build: $(SRCS)
	@echo "Compiling source files..."
	@mkdir -p bin
	$(CC) $(CFLAGS) $^ -o bin/$(EXECUTABLE_NAME) $(shell root-config --libs)
	@echo "Build complete. Executable is at bin/$(EXECUTABLE_NAME)"

clean:
	@echo "Cleaning build artifacts..."
	@if [ -d bin ]; then rm -rf bin; fi
	@echo "Clean complete."