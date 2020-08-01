###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

# directory for placing .d dependency files
DEP_DIR := $(OBJ_DIR)/.deps

##########################################################
DEP_FLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d
## CC COMPILER OPTIONS ##

# CC compiler options:
CC=g++
#-static-libasan -fsanitize=address -fno-omit-frame-pointer -ltcmalloc
CC_FLAGS=-fopenmp -g -O3 -Iinclude -static-libasan -fsanitize=address -fno-omit-frame-pointer
CC_LD_FLAGS=-Wl,--no-as-needed -lprofiler  -Wl,--as-needed
CC_LIBS=

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc -ccbin g++
NVCC_FLAGS=-g -Xcompiler -fopenmp -O3 -Iinclude
NVCC_LD_FLAGS=-lprofiler
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

COMPILE.c = $(CC) $(DEP_FLAGS) $(CC_FLAGS) $(CC_LD_FLAGS) -c
LINK.c = $(CC) $(CC_FLAGS) $(CC_LD_FLAGS)

##########################################################

## Make variables ##

SRCS = $(notdir $(wildcard $(SRC_DIR)/*.cpp))

# Target executable name:
EXE = testMinHash testAligner testContig testConsensus testMyers testAlignAlgs

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
all: $(EXE)

testMinHash : $(addprefix $(OBJ_DIR)/,testMinHash.o NanoporeReads.o)
testAligner : $(addprefix $(OBJ_DIR)/,testAligner.o NanoporeReads.o ReadAligner.o)
testContig : $(addprefix $(OBJ_DIR)/,testContigGenerator.o NanoporeReads.o ReadAligner.o Contig.o)
testConsensus : $(addprefix $(OBJ_DIR)/,testConsensus.o Consensus.o Contig.o NanoporeReads.o ReadAligner.o myers.o Edits.o)
testMyers : $(addprefix $(OBJ_DIR)/,testMyers.o myers.o Edits.o)
testAlignAlgs : $(addprefix $(OBJ_DIR)/,testAlignAlgs.o AlignerTester.o myers.o Edits.o)

$(EXE) :
	$(LINK.c) $^ -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(DEP_DIR)/%.d | $(DEP_DIR)
	$(COMPILE.c) $(OUTPUT_OPTION) $<

$(DEP_DIR): ; @mkdir -p $@

DEP_FILES := $(SRCS:%.cpp=$(DEP_DIR)/%.d)

$(DEP_FILES):
include $(wildcard $(DEP_FILES))

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE) bin/.deps/*

.PHONY: clean all



