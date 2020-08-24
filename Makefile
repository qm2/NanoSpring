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

SUBDIRS := libbsc

##########################################################
DEP_FLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d
## CC COMPILER OPTIONS ##

# CC compiler options: (we need g++ version >= 8)
CC=g++
#-static-libasan -fsanitize=address -fno-omit-frame-pointer -ltcmalloc
CC_FLAGS=-std=c++17 -fopenmp -g -O3 -Iinclude -Ilibbsc \
#-static-libasan -fsanitize=address -fno-omit-frame-pointer
CC_LD_FLAGS=
#-Wl,--no-as-needed -lprofiler
CC_LIBS=-lstdc++fs -Wl,--no-as-needed -lprofiler -Wl,--as-needed

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc -ccbin g++
NVCC_FLAGS=-g -Xcompiler -fopenmp -O3 -Iinclude -Ilibbsc
NVCC_LD_FLAGS=-lprofiler
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

COMPILE.c = $(CC) $(DEP_FLAGS) $(CC_FLAGS) -c
LINK.c = $(CC) $(CC_FLAGS) $(CC_LD_FLAGS)

##########################################################

## Make variables ##

SRCS = $(notdir $(wildcard $(SRC_DIR)/*.cpp))

# Target executable name:
EXE = testMinHash testAligner testContig testCompressor \
	testMyers testAlignAlgs validateResult testDecompressor

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
all: $(SUBDIRS) $(EXE)

testMinHash : $(addprefix $(OBJ_DIR)/,testMinHash.o NanoporeReads.o)
testAligner : $(addprefix $(OBJ_DIR)/,testAligner.o NanoporeReads.o ReadAligner.o)
testContig : $(addprefix $(OBJ_DIR)/,testContigGenerator.o NanoporeReads.o ReadAligner.o Contig.o)
testCompressor : $(addprefix $(OBJ_DIR)/,testCompressor.o Compressor.o Consensus.o Contig.o \
	NanoporeReads.o ReadAligner.o \
	LocalMyersRollBack.o LocalMyersRollBackOld.o LocalMyers.o MyersAligner.o Edits.o\
	bsc.o) \
	$(addprefix libbsc/, libbsc.a)
testMyers : $(addprefix $(OBJ_DIR)/,testMyers.o \
	LocalMyersRollBack.o LocalMyersRollBackOld.o LocalMyers.o MyersAligner.o Edits.o)
testAlignAlgs : $(addprefix $(OBJ_DIR)/,testAlignAlgs.o AlignerTester.o \
	LocalMyersRollBack.o LocalMyersRollBackOld.o LocalMyers.o MyersAligner.o Edits.o)
validateResult : $(addprefix $(OBJ_DIR)/,validateResult.o AlignerTester.o)
testDecompressor : $(addprefix $(OBJ_DIR)/,testDecompressor.o Decompressor.o bsc.o) \
	$(addprefix libbsc/, libbsc.a)
# test : $(addprefix $(OBJ_DIR)/, bsc.o) $(addprefix libbsc/, libbsc.a)

$(EXE) :
	$(LINK.c) $^ -o $@ $(CC_LIBS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(DEP_DIR)/%.d | $(DEP_DIR)
	$(COMPILE.c) $(OUTPUT_OPTION) $<

$(DEP_DIR): ; @mkdir -p $@

DEP_FILES := $(SRCS:%.cpp=$(DEP_DIR)/%.d)

$(DEP_FILES):
include $(wildcard $(DEP_FILES))

$(SUBDIRS):
	$(MAKE) -C $@

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE) bin/.deps/*
	$(MAKE) clean -C $(SUBDIRS)

.PHONY: clean all $(SUBDIRS)



