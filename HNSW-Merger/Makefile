CC = g++ 
CFLAGS_DEBUG = -g
CFLAGS_RELEASE = -O3
CFLAGS_IGNORE = -Wno-write-strings -Wno-unused-result -Wno-format 
CFLAG_FAISS = -I/usr/local/include -L/usr/local/lib -fopenmp
SRC = test_merge.cpp
SRC_AVAIL = availability_test.cpp
TARGET = merge
AVAIL = avail

all: exp

build:
# 	g++ $(CFLAGS_IGNORE) -O3 build_index_on_idx.cpp -o builds $(CFLAG_FAISS)
	g++ $(CFLAGS_IGNORE) -O3 build_index.cpp -o builds $(CFLAG_FAISS)

exp:
	g++ $(CFLAGS_IGNORE) -O3 experiment.cpp -o exps $(CFLAG_FAISS)

debug:
	g++ $(CFLAGS_IGNORE) -g experiment.cpp -o exps $(CFLAG_FAISS)
	gdb ./exps

clean:
	rm -f $(TARGET)
