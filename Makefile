CC ?= clang
BUILD_DIR := build
TARGET := $(BUILD_DIR)/patterson
TEST_TARGET := $(BUILD_DIR)/pointer-equivalence

CPPFLAGS ?=
CFLAGS ?= -Wall -g -O3 -mavx -mavx2 -mtune=native -march=native -ffast-math -funroll-loops
OPENMP_FLAGS ?= -fopenmp
DEPFLAGS := -MMD -MP
LDFLAGS ?=
LDLIBS ?=

SOURCES := main.c Patterson.c ecole.c chash.c lu.c sha3.c inv_mat.c debug.c
OBJECTS := $(SOURCES:%.c=$(BUILD_DIR)/%.o)
TEST_OBJECTS := $(filter-out $(BUILD_DIR)/main.o,$(OBJECTS))

.PHONY: all test clean

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: %.c common.h global.h struct.h simd.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPENMP_FLAGS) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/sha3.o: sha3.h

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OPENMP_FLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/pointer-equivalence.o: tests/pointer_equivalence.c common.h global.h struct.h simd.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPENMP_FLAGS) $(DEPFLAGS) -c $< -o $@

$(TEST_TARGET): $(BUILD_DIR)/pointer-equivalence.o $(TEST_OBJECTS)
	$(CC) $(LDFLAGS) $(OPENMP_FLAGS) $^ $(LDLIBS) -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

-include $(OBJECTS:.o=.d) $(BUILD_DIR)/pointer-equivalence.d

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.d $(TARGET) $(TEST_TARGET)
	rmdir $(BUILD_DIR) 2>/dev/null || true
