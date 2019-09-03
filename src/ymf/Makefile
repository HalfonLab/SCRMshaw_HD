CC = g++ 
CCFLAGS = -Iinclude -I../matrix -g -Wno-deprecated -O2 -Wno-write-strings
OBJ = .obj
SRC = src
BIN = bin
CONFIG = config
PREFIX = ../pac
LIB = $(PREFIX)/lib

install: $(BIN)/statsvar $(BIN)/preproc
#	cp $(BIN)/* $(PREFIX)/bin
#	cp $(CONFIG)/* $(PREFIX)/bin

all: $(BIN)/statsvar $(BIN)/preproc

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CC) $(CCFLAGS) -c $< -o $@  

$(BIN)/statsvar : $(OBJ)/main.o $(OBJ)/expect.o $(OBJ)/variance.o
	$(CXX) $^ -o $@ -lm

$(BIN)/preproc : $(OBJ)/preproc-main.o $(OBJ)/preproc_ymf.o
	$(CXX) -L $(LIB) $^ -o $@ -lmat -lm

clean:
	rm -rf ${OBJ}/*.o ${BIN}/*

