.DEFAULT_GOAL := all


GPP = g++


CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math

ROOT_DIR = checks
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
$(shell mkdir -p $(BIN_DIR))


#DEPS = time_vector_checks 
#FULL_DEPS = $(patsubst %,$(SRC_DIR)/%.cpp,$(DEPS)) #Pad names with dir
##$(info $$OBJ is [${FULL_DEPS}])


all:
	$(GPP) $(CPP_FLAGS) -o $(BIN_DIR)/time_vector_checks $(SRC_DIR)/time_vector_checks.cpp -ljsoncpp
	$(GPP) $(CPP_FLAGS) -o $(BIN_DIR)/initialization_checks $(SRC_DIR)/initialization_checks.cpp -ljsoncpp
	$(GPP) $(CPP_FLAGS) -o $(BIN_DIR)/get_map_path $(SRC_DIR)/get_map_path.cpp -lgerlumph


clean:
	$(RM) -r $(BIN_DIR)/*