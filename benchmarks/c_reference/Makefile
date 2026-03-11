MDTRAJ_XTC := $(HOME)/ghq/github.com/mdtraj/mdtraj/mdtraj/formats/xtc
SRC_DIR := $(MDTRAJ_XTC)/src
INC_DIR := $(MDTRAJ_XTC)/include

CC := cc
CFLAGS := -O3 -I$(INC_DIR)

SRCS := bench.c $(SRC_DIR)/xdrfile.c $(SRC_DIR)/xdrfile_xtc.c $(SRC_DIR)/xdrfile_trr.c $(SRC_DIR)/xdr_seek.c
TARGET := bench

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

run: $(TARGET)
	cd ../.. && benchmarks/c_reference/bench

clean:
	rm -f $(TARGET)
