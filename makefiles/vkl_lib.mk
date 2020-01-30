GPP = g++
FORTRAN = gfortran


############## START: VKL_LIB
FORTRAN_FLAGS = -fPIC -frounding-math
CPP_FLAGS = -std=c++11 -fPIC -g -frounding-math
CPP_LIBS  = -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL

ROOT_DIR = common/vkl_lib
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/inc
LIB_DIR = $(ROOT_DIR)/lib
OBJ_DIR = $(ROOT_DIR)/obj
$(shell mkdir -p $(OBJ_DIR))
$(shell mkdir -p $(LIB_DIR))


DEPS = imagePlane.hpp massModels.hpp sourceProfile.hpp sourcePlane.hpp nonLinearPars.hpp tableDefinition.hpp covKernels.hpp
OBJ  = imagePlane.o   massModels.o   sourceProfile.o   sourcePlane.o   nonLinearPars.o   fastell.o
FULL_DEPS = $(patsubst %,$(INC_DIR)/%,$(DEPS)) #Pad names with dir
FULL_OBJ  = $(patsubst %,$(OBJ_DIR)/%,$(OBJ))  #Pad names with dir
#$(info $$OBJ is [${FPROJECT_DEPS}])

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(FULL_DEPS)
	$(GPP) $(CPP_FLAGS) -I $(INC_DIR) -c -o $@ $<
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FORTRAN) $(FORTRAN_FLAGS) -c -o $@ $<

vkl_lib: $(FULL_OBJ)
	$(GPP) -shared -Wl,-soname,libvkl.so -o $(LIB_DIR)/libvkl.so $(FULL_OBJ) $(CPP_LIBS)
vkl_lib_clean:
	$(RM) -r $(OBJ_DIR)/* $(LIB_DIR)/*
############## END: VKL_LIB
