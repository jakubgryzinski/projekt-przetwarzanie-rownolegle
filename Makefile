NVCC = nvcc
TARGET = main
SRC_DIR = src
BIN_DIR = bin
OUT_DIR = outputs
ARCH = -gencode arch=compute_61,code=sm_61
ARGS =

$(shell mkdir -p $(BIN_DIR) $(OUT_DIR))

all: $(TARGET)

$(TARGET): $(SRC_DIR)/main.cu
	@$(NVCC) $(ARCH) --std=c++11 -o $(BIN_DIR)/$(TARGET) $(SRC_DIR)/main.cu

run: $(TARGET)
	@./$(BIN_DIR)/$(TARGET) $(ARGS)

test: $(SRC_DIR)/main.cpp
	@g++ --std=c++11 -o $(BIN_DIR)/test $(SRC_DIR)/main.cpp
	@./$(BIN_DIR)/test $(ARGS)

debug: $(SRC_DIR)/main.cu
	@$(NVCC) $(ARCH) --std=c++11 -g -G -o $(BIN_DIR)/$(TARGET) $(SRC_DIR)/main.cu

clean:
	@rm -f $(BIN_DIR)/$(TARGET) $(BIN_DIR)/test *.o
	@rm -f $(OUT_DIR)/*

.PHONY: all run test debug clean