CXX = g++
CXXFLAGS = -O -w -Iinclude
OBJ = .obj
SRC = src
BIN = bin
PREFIX = ../..

install: $(BIN)/compare_hex $(BIN)/word_score
	cp $(BIN)/* $(PREFIX)/bin

all: $(BIN)/compare_hex $(BIN)/word_score 

$(OBJ)/%.o: $(SRC)/%.cpp
#	$(CXX) $(CXXFLAGS) $(SQLCXXFLAGS) -c $< -o $@ 
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/compare_hex : $(OBJ)/compare_hex.o $(OBJ)/listop.o
	$(CXX) $^ -o $@

$(BIN)/word_score : $(OBJ)/word_score.o $(OBJ)/listop.o
	$(CXX) $^ -o $@

clean : 
	rm -rf ${OBJ}/*.o ${BIN}/*
