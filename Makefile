###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda

##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=g++
CC_FLAGS=-fopenmp -g -O3
CC_LD_FLAGS=
CC_LIBS=

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc
NVCC_FLAGS=-g -Xcompiler -fopenmp -O3
NVCC_LD_FLAGS=-lprofiler
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

##########################################################

## Make variables ##

# Target executable name:
EXE = testMinHash testAligner

# Object files:
OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/cuda_kernel.o

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
all: $(EXE)

testMinHash : $(OBJ_DIR)/testMinHash.o $(OBJ_DIR)/NanoporeReads.o
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LD_FLAGS) $^ -o $@

testAligner : $(OBJ_DIR)/testAligner.o $(OBJ_DIR)/NanoporeReads.o $(OBJ_DIR)/ReadAligner.o
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LD_FLAGS) $^ -o $@

# Compile main .cpp file to object files:
$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp include/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile CUDA source files to object files:

$(OBJ_DIR)/testAligner.o: $(SRC_DIR)/testAligner.cu $(INC_DIR)/testAligner.cuh $(INC_DIR)/NanoporeReads.cuh $(INC_DIR)/ReadAligner.cuh
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(OBJ_DIR)/testMinHash.o: $(SRC_DIR)/testMinHash.cu $(INC_DIR)/testMinHash.cuh $(INC_DIR)/NanoporeReads.cuh
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(OBJ_DIR)/NanoporeReads.o : $(SRC_DIR)/NanoporeReads.cu $(INC_DIR)/NanoporeReads.cuh
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(OBJ_DIR)/ReadAligner.o : $(SRC_DIR)/ReadAligner.cu $(INC_DIR)/ReadAligner.cuh $(INC_DIR)/NanoporeReads.cuh
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)

.PHONY: clean all



