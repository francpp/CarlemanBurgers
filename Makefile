# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Iinclude -Iinclude/nlohmann -Isrc -I${mkEigenInc}

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp) \
            $(wildcard $(SRC_DIR)/params/*.cpp) \
            $(wildcard $(SRC_DIR)/discretization/*.cpp) \
            $(wildcard $(SRC_DIR)/initial_conditions/*.cpp) \
            $(wildcard $(SRC_DIR)/matrix/*.cpp) \
            $(wildcard $(SRC_DIR)/solvers/*.cpp) \
            $(wildcard $(SRC_DIR)/error_analysis/*.cpp) \
            $(wildcard $(SRC_DIR)/utils/*.cpp) \
            main.cpp
OBJ_FILES = $(SRC_FILES:%.cpp=$(BUILD_DIR)/%.o)

# Target executable
TARGET = main_executable

# Default rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $(TARGET)

# Compilation
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Phony targets
.PHONY: all clean
