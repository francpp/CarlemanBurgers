# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Isrc

# Directories
SRC_DIR = src
BUILD_DIR = build

# Optional include directories (conditionally added)
ifneq (,$(shell grep -rl '#include "json.hpp"' $(SRC_DIR)))
CXXFLAGS += -Iinclude/nlohmann
endif

ifneq (,$(shell grep -rl '#include <Eigen/' $(SRC_DIR)))
CXXFLAGS += -I${mkEigenInc}
endif

ifneq (,$(shell grep -rl '#include "GetPot"' $(SRC_DIR)))
CXXFLAGS += -Iinclude
endif

ifneq (,$(shell grep -rl '#include "boost/' $(SRC_DIR)))
CXXFLAGS += -I${mkBoostInc}
LDFLAGS += -L${mkBoostLib} -lboost_iostreams -lboost_system -lboost_filesystem
endif

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp) \
            $(wildcard $(SRC_DIR)/*/*.cpp) \
            main.cpp

OBJ_FILES = $(SRC_FILES:%.cpp=$(BUILD_DIR)/%.o)

# Target executable
TARGET = main_executable

# Default rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $(TARGET) $(LDFLAGS)

# Compilation
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Phony targets
.PHONY: all clean
