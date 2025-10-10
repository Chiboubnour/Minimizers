SHELL := /bin/bash

# Configuration
TARGET ?= hw
PLATFORM ?= xilinx_u280_gen3x16_xdma_1_202211_1
KERNEL_NAME := krnl_hash_simple

# Directories
BUILD_DIR := build_dir.$(TARGET).$(PLATFORM)
SRC_DIR := src

# Files
HOST_SRC := $(SRC_DIR)/host_simple.cpp
KERNEL_SRC := $(SRC_DIR)/$(KERNEL_NAME).cpp
CONFIG_FILE := $(KERNEL_NAME).cfg

# Output files
XO_FILE := $(BUILD_DIR)/$(KERNEL_NAME).xo
XCLBIN_FILE := $(BUILD_DIR)/$(KERNEL_NAME).xclbin
HOST_EXE := host_simple

# Tools
VPP := v++
CXX := g++

# Compiler flags
CXXFLAGS := -std=c++17 -O2 -Wall -g -I$(XILINX_XRT)/include
LDFLAGS := -L$(XILINX_XRT)/lib -lxrt_coreutil -luuid -lpthread

# Vitis flags
VPP_FLAGS := -t $(TARGET) --platform $(PLATFORM) --save-temps -g

# Default target
all: build

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile kernel to XO
$(XO_FILE): $(KERNEL_SRC) $(CONFIG_FILE) | $(BUILD_DIR)
	$(VPP) -c $(VPP_FLAGS) \
		-k $(KERNEL_NAME) \
		--temp_dir $(BUILD_DIR) \
		--output $(XO_FILE) $<

# Link XO to XCLBIN
$(XCLBIN_FILE): $(XO_FILE)
	$(VPP) -l $(VPP_FLAGS) \
		--config $(CONFIG_FILE) \
		--temp_dir $(BUILD_DIR) \
		--output $(XCLBIN_FILE) $(XO_FILE)

# Compile host
$(HOST_EXE): $(HOST_SRC)
	$(CXX) $(CXXFLAGS) -o $(HOST_EXE) $(HOST_SRC) $(LDFLAGS)

# Build everything
build: $(XCLBIN_FILE) $(HOST_EXE)

# Run on board
run: build
	./$(HOST_EXE) $(XCLBIN_FILE) 0 512

# Run with custom data size
run-custom: build
	./$(HOST_EXE) $(XCLBIN_FILE) 0 $(DATA_SIZE)

# Clean generated files
clean:
	rm -rf $(BUILD_DIR) $(HOST_EXE) *.log *.info *.wdb

# Full cleanup
cleanall: clean
	rm -rf _x* .Xil package.* *.csv *.json *.protoinst

# Help
help:
	@echo "=== Makefile pour Hash Simple ==="
	@echo ""
	@echo "Targets disponibles:"
	@echo "  all        - Build tout (défaut)"
	@echo "  build      - Build kernel et host"
	@echo "  run        - Build et run avec 512MB de données"
	@echo "  run-custom - Build et run avec taille personnalisée (DATA_SIZE=xxx)"
	@echo "  clean      - Nettoyer les fichiers générés"
	@echo "  cleanall   - Nettoyage complet"
	@echo ""
	@echo "Variables:"
	@echo "  TARGET     - Cible de build (hw, hw_emu) [défaut: hw]"
	@echo "  PLATFORM   - Nom de la plateforme [défaut: xilinx_u280_gen3x16_xdma_1_202211_1]"
	@echo "  DATA_SIZE  - Taille des données en MB pour run-custom"
	@echo ""
	@echo "Exemples:"
	@echo "  make run TARGET=hw_emu"
	@echo "  make run-custom DATA_SIZE=256"

.PHONY: all build run run-custom clean cleanall help
