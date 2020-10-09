#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility> 
using namespace std;

//Create matrix for optimal alignment, write it
vector<vector<int>> scoreMatrix(string seq1, string seq2);
void printOptAlignment(vector<vector<int>>& matrix, string seq1, string seq2);
//Write optimal score
//Check for multiple optimal alignments
//Find best optimal alignment and write the best alignment out on two lines
//EC: Write out all possible optimal alignments

/*SCORING FUNCTION:
2 for match, -1 for mismatch, -2 match against gap
***************************************************/

const int MATCH = 2;
const int MISMATCH = -1;
const int GAP = -2;


int main(int argc, char** argv) {
  string inLine;
  string seq1, seq2;
  bool firstSeq = true;

  cout << "Enter File Name" << endl;
  getline(cin, inLine); // Get file name
  // Open file
  ifstream infile(inLine);
  if (infile.is_open())
  {
    while ( getline (infile,inLine) )
    {
      if ( inLine.size() && inLine.at(0) == '>') {
        getline( infile, (firstSeq ? seq1 : seq2));
        firstSeq = !firstSeq;
      }
    }
    infile.close();
  }
  else cout << "Unable to open file" << endl;
  
  cout << "Comparing " << seq1 << " and " << seq2 << endl << endl;

  auto matrix = scoreMatrix(seq1, seq2);
  
  ofstream outfile;
  // Output the optimal score
  outfile.open("assignment1.o1");
  if (outfile.is_open())
  {
    cout << "Optimal sequence alignment score is: " << matrix.back().back() << endl << endl;
    outfile << matrix.back().back() << endl;
    outfile.close();
  }
  else {
    cout << "Unable to open file";
    return -1;
  }

  // Output the score matrix
  outfile.open("assignment1.o2");
  if (outfile.is_open())
  {
    cout << "Score Matrix:" << endl;
    for (auto row : matrix) {
      for (auto num : row) {
        outfile << num << ' ';
        cout << num << ' ';
      }
      outfile << endl;
      cout << endl;
    }
    outfile.close();
  }
  else cout << "Unable to open file";
  printOptAlignment(matrix, seq1, seq2);
  // Set the two sequences
  return 0;
}

vector<vector<int>> scoreMatrix(string seq1, string seq2) {
  //vector<vector<int>> matrix(m.size() + 1, n.size() + 1);
  //vector<pair<int,int>> path;
  vector<vector<int>> matrix(seq1.size() + 1, vector<int> (seq2.size()+1,0));
  for (int i = 0; i < seq1.size() + 1; i++) {
      matrix[i][0] = i * GAP;
  }

  for (int j = 0; j < seq2.size() + 1; j++) {
    matrix[0][j] = j * GAP;
  }

  for (int i = 1; i < seq1.size() + 1; i++) {
    for (int j = 1; j < seq2.size() + 1; j++) {      
      
      if (seq1[i-1] == seq2[j-1]) {
        matrix[i][j] = max({matrix[i-1][j-1] + MATCH,
                            matrix[i-1][j] + GAP, 
                            matrix[i][j-1] + GAP});
      }
      else { // seq1[i-1] != seq2[j-1]
        matrix[i][j] = max({matrix[i-1][j-1] + MISMATCH,
                            matrix[i-1][j] + GAP,
                            matrix[i][j-1] + GAP});
      }
      
    }
  }
  
  return matrix;
}

void printOptAlignment(vector<vector<int>>& matrix, string seq1, string seq2){
  
  int i = matrix.size() - 1;
  int j = matrix[0].size() - 1;
  //I am currently testing things with print statements
  while (i > 0 && j > 0) { //OR breaks..    
    if (matrix[i-1][j-1] == matrix[i][j] - MATCH
      || matrix[i-1][j-1] == matrix[i][j] - MISMATCH) {
        //cout << i <<" and " << j <<"\n";
        i--;
        j--;        
    }    
    //cout << "after " << i <<" and " << j <<"\n";
    if (i > 0 && matrix[i-1][j] == matrix[i][j] - GAP) {
      seq2.insert(i-1, "-");
      //cout << "gap " << i <<" and " << j <<"\n";
      //cout << matrix[i-1][j] << " - 2 = " << matrix[i][j] << "\n";
      i--;
      //cout << "after gap " << i <<" and " << j <<"\n";
    }
    if (j > 0 && matrix[i][j-1] == matrix[i][j] - GAP) {
      seq1.insert(j-1,"-");
      j--;
    }
  }
  cout << seq1 << endl;
  cout << seq2 << endl;
}



// opt(matrix, seq1, seq2) {
  

// }