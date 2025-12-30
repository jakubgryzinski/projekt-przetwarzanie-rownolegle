NVCC = nvcc
TARGET = main
ARCH = -gencode arch=compute_61,code=sm_61

all: $(TARGET)

$(TARGET): main.cu
	$(NVCC) $(ARCH) --std=c++11 -o $(TARGET) main.cu

run: $(TARGET)
	./$(TARGET)

debug: main.cu
	$(NVCC) $(ARCH) --std=c++11 -g -G -o $(TARGET) main.cu

clean:
	rm -f $(TARGET) *.o

.PHONY: all run debug clean