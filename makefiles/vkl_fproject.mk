GPP = g++
MPICC = mpic++
FORTRAN = gfortran


############## START: Fproject
FORTRAN_FLAGS = -fPIC -frounding-math
FPROJECT_FLAGS = -std=c++11 -fPIC -g -frounding-math
FPROJECT_LIBS  = -lgfortran -lCCfits -lcfitsio -ljsoncpp -lgmp -lCGAL

FPROJECT_DIR = lensed_extended_source/vkl_fproject
SRC_DIR = $(FPROJECT_DIR)/src
INC_DIR = $(FPROJECT_DIR)/inc
BIN_DIR = $(FPROJECT_DIR)/bin
OBJ_DIR = $(FPROJECT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(BIN_DIR))


FP_DEPS = initFuncs.hpp polyclip.hpp polygons.hpp sourceProfile.hpp massModels.hpp sourcePlane.hpp imagePlane.hpp nonLinearPars.hpp contour-tracing.hpp tableDefinition.hpp covKernels.hpp
FP_OBJ  = initFuncs.o   polyclip.o   polygons.o   sourceProfile.o   massModels.o   sourcePlane.o   imagePlane.o   nonLinearPars.o   contour-tracing.o   fastell.o           fproject.o
FPROJECT_DEPS = $(patsubst %,$(INC_DIR)/%,$(FP_DEPS)) #Pad names with dir
FPROJECT_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(FP_OBJ))  #Pad names with dir
#$(info $$OBJ is [${FPROJECT_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FPROJECT_DEPS)
	$(GPP) $(FPROJECT_FLAGS) -I $(INC_DIR) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FORTRAN) $(FORTRAN_FLAGS) -c -o $@ $<

fproject: $(FPROJECT_OBJ)
	$(GPP) $(FPROJECT_FLAGS) -I $(INC_DIR) -o $(BIN_DIR)/fproject $(FPROJECT_OBJ) $(FPROJECT_LIBS)
fproject_clean:
	$(RM) -r $(OBJ_DIR)/* $(BIN_DIR)/*
############## END: Fproject
