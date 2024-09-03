# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -Isrc

# Directories
SRC_DIR = src
BUILD_DIR = build
DOCS_DIR = docs
INCLUDE_DIR = include

# PACS directories
PACS_ROOT = /home/fpettenon/PACS/pacs-examples
PACS_INC_DIR = $(PACS_ROOT)/Examples/src/Utilities/
PACS_EXTRAS_DIR = $(PACS_ROOT)/Extras/
PACS_MUPARSER_DIR = $(PACS_EXTRAS_DIR)/muparser/include/
PACS_JSON_DIR = $(PACS_EXTRAS_DIR)/json/include/nlohmann/

# Optional include directories (conditionally added)
ifneq (,$(shell grep -rl '#include "json.hpp"' $(SRC_DIR)))
CXXFLAGS += -I$(INCLUDE_DIR)/nlohmann
endif

ifneq (,$(shell grep -rl '#include <Eigen/' $(SRC_DIR)))
CXXFLAGS += -I${mkEigenInc}
endif

ifneq (,$(shell grep -rl '#include "GetPot"' $(SRC_DIR)))
CXXFLAGS += -I$(INCLUDE_DIR)
endif

ifneq (,$(shell grep -rl '#include "boost/' $(SRC_DIR)))
CXXFLAGS += -I${mkBoostInc}
LDFLAGS += -L${mkBoostLib} -lboost_iostreams -lboost_system -lboost_filesystem
endif

ifneq (,$(shell grep -rl '"utils/muparser_fun.hpp"' $(SRC_DIR)))
CXXFLAGS += -I$(INCLUDE_DIR)
LDFLAGS += -L/home/fpettenon/PACS/pacs-examples/Examples/lib -lmuparser
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

# Generate documentation
docs:
	doxygen Doxyfile

# Clean rule
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(DOCS_DIR)
	find $(INCLUDE_DIR) -type f ! -name '.gitignore' -delete
	find $(INCLUDE_DIR) -type d -empty -delete

# Install headers and other files
install:
	# Copy GetPot and gnuplot headers
	cp $(PACS_ROOT)/Examples/src/Utilities/GetPot $(INCLUDE_DIR)/
	cp $(PACS_ROOT)/Examples/src/Utilities/gnuplot-iostream.hpp $(INCLUDE_DIR)/
	# Copy json headers
	cp -r $(PACS_JSON_DIR) $(INCLUDE_DIR)/
	# Copy muparser headers
	cp $(PACS_MUPARSER_DIR)*.h $(INCLUDE_DIR)/

# Phony targets
.PHONY: all clean docs install
