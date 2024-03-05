# Builds Mini Chem executable

# Application name - oasis
target_exec := mini-che
# Build application folder
build_dir := ./build
# Source code folder
src_dir := ./src

#
cuda := /usr/local/cuda

# Streaming Multiprocessor version
gencode_sm35	:= -gencode arch=compute_35,code=sm_35 # Kepler
gencode_sm37	:= -gencode arch=compute_37,code=sm_37 # Kepler
gencode_sm50	:= -gencode arch=compute_50,code=sm_50 # Maxwell
gencode_sm52	:= -gencode arch=compute_52,code=sm_52 # Maxwell
gencode_sm53	:= -gencode arch=compute_53,code=sm_53 # Maxwell
gencode_sm60	:= -gencode arch=compute_60,code=sm_60 # Pascal
gencode_sm61	:= -gencode arch=compute_61,code=sm_61 # Pascal
gencode_sm70	:= -gencode arch=compute_70,code=sm_70 # Volta
gencode_sm75  := -gencode arch=compute_75,code=sm_75 # Turing
gencode_sm80  := -gencode arch=compute_80,code=sm_80 # Ampere A100
gencode_sm86  := -gencode arch=compute_86,code=sm_86 # Ampere RTX 30xx
gencode_flags := $(gencode_sm86)

# Compiler
cudac:= nvcc

# Find all the CU files we want to compile
src := $(shell find $(src_dir) -name '*.cu')
objs := $(src:%=$(build_dir)/%.o)
deps := $(objs:.o=.d)

src.cu := $(filter $(src_dir)/%.cu, $(src))
objs.cu := $(src.cu:%=$(build_dir)/%.o)

# Every folder passed to nvcc so that it can find header files
inc_dirs := $(shell find $(src_dir) -type d)
inc_flags := $(addprefix -I,$(inc_dirs))

cuflags := $(inc_flags) $(gencode_flags) -Wno-deprecated-gpu-targets

# The final build step
$(target_exec): $(objs.cu)
	$(cudac) $(cuflags) -o $@ $(objs.cu)

# Build step for CU source
$(build_dir)/%.cu.o: %.cu
	mkdir -p $(dir $@)
	$(cudac) $(cuflags) -dc $< -o $@

.PHONY: clean
clean:
	rm -r $(build_dir)
	rm -f $(target_exec)

# Include the .d makefiles
-include $(deps)
