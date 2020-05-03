.DEFAULT_GOAL := match_to_gerlumph

GPP = g++
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -lsqlite3 -ljsoncpp

ROOT_DIR = variability/extrinsic/match_to_gerlumph
INC_DIR = $(ROOT_DIR)/inc
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))

DEPS = sql_callback.hpp
OBJ  = match_to_gerlumph.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FULL_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -c -o $@ $<

match_to_gerlumph: $(FULL_OBJ)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -o $(BIN_DIR)/match_to_gerlumph $(FULL_OBJ) $(CPP_LIBS)
clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
