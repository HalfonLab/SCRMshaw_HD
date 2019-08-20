CXX = g++
CXXFLAGS = -g -Wno-deprecated -O2 -Wno-write-strings -w 
OBJ = .obj
SRCHEX = ./src/hexmcd/src
SRCIMM = ./src/imm/src
SRCPAC = ./src/pac
SRCYMF = ./src/ymf/src
SRCMATRIX = ./src/matrix
BIN = ./bin
CODE = ./code
LIB = ./src/pac/lib
LIBS = -L$(LIB) -lmat -lm
OBJS = $(OBJ)/stats.o $(OBJ)/preproc.o $(OBJ)/listop.o

install: $(BIN)/compare_hex $(BIN)/word_score $(BIN)/imm_build $(BIN)/imm_score $(BIN)/pac $(BIN)/statsvar $(BIN)/preproc

all: $(BIN)/compare_hex $(BIN)/word_score $(BIN)/imm_build $(BIN)/imm_score $(BIN)/pac $(BIN)/statsvar $(BIN)/preproc

$(OBJ)/%.o: $(SRCHEX)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I ./src/hexmcd/include

$(OBJ)/%.o: $(SRCIMM)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I ./src/imm/include

$(OBJ)/%.o: $(SRCPAC)/src/%.cpp
	$(CXX) $(CXXFLAGS) -D_ADDITIVE -c $< -o $@ -I ./src/pac/include

$(OBJ)/%.o: $(SRCYMF)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I ./src/ymf/include -I ./src/matrix

$(OBJ)/%.o: $(SRCMATRIX)/%.cpp
	$(CXX) -c $< -o $@

$(BIN)/compare_hex : $(OBJ)/compare_hex.o $(OBJ)/listop.o
	$(CXX) $^ -o $@

$(BIN)/word_score : $(OBJ)/word_score.o $(OBJ)/listop.o
	$(CXX) $^ -o $@

$(BIN)/imm_score : $(OBJ)/imm_score.o ${OBJ}/sequence.o ${OBJ}/tree_node.o ${OBJ}/tree.o
	$(CXX) $^ -o $@

$(BIN)/imm_build : $(OBJ)/imm_build.o ${OBJ}/sequence.o ${OBJ}/tree_node.o ${OBJ}/tree.o
	$(CXX) $^ -o $@

$(LIB)/libmat.a : $(OBJ)/matrix.o
	ar rv $@ $^

$(BIN)/statsvar : $(OBJ)/main.o $(OBJ)/expect.o $(OBJ)/variance.o
	$(CXX) $^ -o $@ -lm

$(BIN)/preproc : $(OBJ)/preproc-main.o $(OBJ)/preproc_ymf.o
	$(CXX) -L $(LIB) $^ -o $@ -lmat -lm

$(SRCMATRIX)/libmat.a : $(SRCMATRIX)/matrix.o
	ar rv $(SRCMATRIX)/libmat.a $(SRCMATRIX)/matrix.o

$(BIN)/pac : $(OBJ)/compare_poisson_subset.o $(OBJS) $(LIB)/libmat.a
	$(CXX) -D_ADDITIVE $(OBJ)/compare_poisson_subset.o $(OBJS) -o $@ $(LIBS)

clean : 
	rm -rf ${OBJ}/*.o ${BIN}/* $(LIB)/*
