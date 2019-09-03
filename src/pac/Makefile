CC = g++ 
CCFLAGS = -Iinclude -g -Wno-deprecated -Wno-write-strings
OBJ = .obj
SRC = src
BIN = bin
LIB = lib
LIBS = -L$(LIB) -lmat -lm 
OBJS = $(OBJ)/stats.o $(OBJ)/preproc.o $(OBJ)/listop.o
PREFIX = ../..

install: $(BIN)/pac
	cp $(BIN)/* $(PREFIX)/bin

all: $(BIN)/pac

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CC) $(CCFLAGS) -D_ADDITIVE -c $< -o $@  

$(LIB)/libmat.a : $(OBJ)/matrix.o
	ar rv $@ $^

$(BIN)/pac : $(OBJ)/compare_poisson_subset.o $(OBJS) $(LIB)/libmat.a
	$(CC) -D_ADDITIVE $(OBJ)/compare_poisson_subset.o $(OBJS) -o $@ $(LIBS)

clean:
	rm -rf ${OBJ}/*.o ${BIN}/* ${LIB}/*
