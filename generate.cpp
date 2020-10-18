#include <iostream>
#include <fstream>

using namespace std;

char randChar() {
  char chars[4] = {'A', 'T', 'C', 'G'};
  return chars[rand() % 4];
}

int main(int argc, char** argv) {
  int length;
  string inLine;
  cout << "Enter File Name" << endl;
	getline(cin, inLine); // Get file name
	// Open file
	ofstream outfile(inLine);
	if (!outfile.is_open()) {
    cout << "Could not open file" << endl;
    return -1;
  }

  cout << "Enter length for each sequence" << endl;
  cin >> length;

  outfile << "> Sequence 1" << endl;
  for (int i = 0; i < length; i++) outfile << randChar();
  outfile << endl;
  outfile << endl;
  outfile << "> Sequence 2" << endl;
  for (int i = 0; i < length; i++) outfile << randChar();
  outfile << endl;
  outfile.close();

  cout << "File " << inLine << " has been written to with sequences of length " << length << endl;
}