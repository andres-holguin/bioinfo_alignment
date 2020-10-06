#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

//Create matrix for optimal alignment, write it
vector<vector<int>> scoreMatrix(string m, string n);
//Write optimal score
//Check for multiple optimal alignments
//Find best optimal alignment and write the best alignment out on two lines
//EC: Write out all possible optimal alignments

/*SCORING FUNCTION:
2 for match, -1 for mismatch, -2 match against gap
***************************************************/
struct score {
  const int MATCH = 2;
  const int MISMATCH = -1;
  const int GAP = -2;
};

int main(int argc, char** argv) {
  string line;

  // Open file
  ifstream myfile("sample.fasta");
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      if ( line.size() && line.at(0) == '>') {
        cout << line << endl;
      }
    }
    myfile.close();
  }
  else cout << "Unable to open file" << endl;

  auto matrix = scoreMatrix("hello", "world");
  for (auto row : matrix) {
    for (auto num : row) {
      cout << num << ' ';
    }
    cout << endl;
  }

  // Set the two sequences
  return 0;
}

vector<vector<int>> scoreMatrix(string m, string n) {
  //vector<vector<int>> matrix(m.size() + 1, n.size() + 1);
  vector<vector<int>> matrix(m.size() + 1, vector<int> (n.size()+1,0));
  for (int i = 0; i < n.size() + 1; i++) {
      matrix[i][0] = i * -2;
  }
  

  for (int j = 0; j < m.size()+1; j++) {
    matrix[0][j] = j * -2;
  }
  
  return matrix;
}