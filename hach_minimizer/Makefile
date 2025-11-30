CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -march=native -std=c++17 -I./include -I./src -I./hls
LDFLAGS = -fopenmp
BIN_DIR = ./bin
SRC_DIR = ./src
TEST_DIR = ./test
HLS_DIR = ./hls
DATASET_DIR = ./dataset

# Noms des exécutables
TARGETS = $(BIN_DIR)/minimizer_mono $(BIN_DIR)/minimizer_multi $(BIN_DIR)/minimizer_hls $(BIN_DIR)/generate_minimizers

all: $(BIN_DIR) $(TARGETS)

# Créer le répertoire bin s'il n'existe pas
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Règle pour la version mono-thread
$(BIN_DIR)/minimizer_mono: $(TEST_DIR)/main_mono.o $(SRC_DIR)/hach_minimizer_mono.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Règle pour la version multi-thread
$(BIN_DIR)/minimizer_multi: $(TEST_DIR)/main_multi.o $(SRC_DIR)/hach_minimizer_multi.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Règle pour la version HLS
$(BIN_DIR)/minimizer_hls: $(TEST_DIR)/main_hls.o $(HLS_DIR)/hach_minimizer_hls.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Règle pour le générateur de données
$(BIN_DIR)/generate_minimizers: $(DATASET_DIR)/main_data.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Règles de compilation des fichiers objets
$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(HLS_DIR)/%.o: $(HLS_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DATASET_DIR)/%.o: $(DATASET_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Nettoyage
clean:
	rm -f $(BIN_DIR)/*
	rm -f $(TEST_DIR)/*.o
	rm -f $(SRC_DIR)/*.o
	rm -f $(HLS_DIR)/*.o
	rm -f $(DATASET_DIR)/*.o

.PHONY: all clean