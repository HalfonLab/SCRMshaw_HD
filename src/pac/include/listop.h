#include <list>
#include <vector>
#include <string>

using std::string;

// Split doesnt create memory for the tokens returned,
// only creates memory for the list holding them.
std::list<char *> *Split(char *toks, char * str);

// Join creates memory for the string returned
char *Join(char *tok, std::list<char *> *l);
char *Join(char *tok, std::list<string> *l);

// Shift doesnt delete the memory taken by the head object 
char *Shift(std::list<char *> *l);

// These will create space for the list returned
std::list<int> *StringListToIntList(std::list<char *> *l);
std::vector<int> *StringListToIntVector(std::list<char *> *l);

// regular util functions
char* itoa( int value, char* result, int base );
void Chomp(char *line);
