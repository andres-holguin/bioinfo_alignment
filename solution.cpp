#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <iomanip>
using namespace std;

vector<vector<int>> scoreMatrix(string seq1, string seq2);
void printIfMultiple(bool ifMulti);
void printOptAlignment(vector<vector<int>>& matrix, vector<vector<int*>>& path, string seq1, string seq2);
//EC: Write out all possible optimal alignments

/*************SCORING FUNCTION:**********************/
/*2 for match, -1 for mismatch, -2 match against gap*/
/****************************************************/
const int MATCH = 2;
const int MISMATCH = -1;
const int GAP = -2;

int main(int argc, char** argv) {
	string inLine;
	string seq1, seq2;
	bool firstSeq = true;

	//To read in argument, if possible: inLine = argv[1];

	cout << "Enter File Name" << endl;
	getline(cin, inLine); // Get file name
	// Open file
	ifstream infile(inLine);
	if (infile.is_open())
	{
		while (getline(infile, inLine))
		{
			if (inLine.size() && inLine.at(0) == '>') {
				getline(infile, (firstSeq ? seq1 : seq2));
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
		cout << "Unable to open file.\n";
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
				cout << right << setw(4) << num << ' ';
			}
			outfile << endl;
			cout << endl;
		}
		outfile.close();
	}
	else cout << "Unable to open file.\n";

	return 0;
}

//Sets up the score matrix using the recurrence equation
vector<vector<int>> scoreMatrix(string seq1, string seq2) {
  //Matrix, path, and score variables
	vector<vector<int>> matrix(seq1.size() + 1, vector<int>(seq2.size() + 1, 0));
  vector<vector<int*>> path(seq1.size() + 1, vector<int*>(seq2.size() + 1));
	int match, mismatch, gap2, gap1;

  //Initial setup
	matrix[0][0] = 0;
	path[0][0] = &matrix[0][0]; //To avoid potential errors, the first element will point to itself
  
	for (int i = 1; i < seq1.size() + 1; i++) {
		matrix[i][0] = i * GAP;
		path[i][0] = &matrix[i - 1][0];
	}
	for (int j = 1; j < seq2.size() + 1; j++) {
		matrix[0][j] = j * GAP;
		path[0][j] = &matrix[0][j - 1];
	}
	//The rest of the matrix
	for (int i = 1; i < seq1.size() + 1; i++) {
		for (int j = 1; j < seq2.size() + 1; j++) {
			match = matrix[i - 1][j - 1] + MATCH;
			mismatch = matrix[i - 1][j - 1] + MISMATCH;
			gap2 = matrix[i - 1][j] + GAP;
			gap1 = matrix[i][j - 1] + GAP;

			if (seq1[i - 1] == seq2[j - 1]) {
				matrix[i][j] = max({ match, gap2, gap1 });
			}
			else { // seq1[i-1] != seq2[j-1]
				matrix[i][j] = max({ mismatch, gap2, gap1 });
			}

			//Storing info for path
			if (matrix[i][j] == match || matrix[i][j] == mismatch) { //Diagonal
				path[i][j] = &matrix[i - 1][j - 1];				
			}
			if (matrix[i][j] == gap2) { //Up
				path[i][j] = &matrix[i - 1][j];				
			}
			if (matrix[i][j] == gap1) { //Left
				path[i][j] = &matrix[i][j - 1];				
			}			
		}
	}	
	printOptAlignment(matrix, path, seq1, seq2);
	return matrix;
}

//Goes through the path and prints an optimal alignment
void printOptAlignment(vector<vector<int>>& matrix, vector<vector<int*>>& path, string seq1, string seq2) {
	ofstream alignFile;
	string newSeq1 = "";
	string newSeq2 = "";
	int i = matrix.size() - 1;
	int j = matrix[0].size() - 1;
	int match, mismatch, gap2, gap1;
	bool multiplePaths = false;
	
	while (i > 0 && j > 0) {
		int multi = 0;
		match = matrix[i - 1][j - 1] + MATCH;
		mismatch = matrix[i - 1][j - 1] + MISMATCH;
		gap2 = matrix[i - 1][j] + GAP;
		gap1 = matrix[i][j - 1] + GAP;

    //Moving diagonally
		if (path[i][j] == &matrix[i - 1][j - 1]) {
			newSeq1.insert(0, 1, seq1[i - 1]);
			newSeq2.insert(0, 1, seq2[j - 1]);
			i--;
			j--;
      //Checking for additional paths
			if (matrix[i][j] == gap2) //Also up
				multi++;
			if (matrix[i][j] == gap1)  //Also left
				multi++;			
		}
    //Moving left
		if (j > 0 && path[i][j] == &matrix[i][j - 1]) {
			newSeq1.insert(0, 1, '-');
			newSeq2.insert(0, 1, seq2[j - 1]);
			j--;
      //Checking for additional paths
			if (matrix[i][j] == match || matrix[i][j] == mismatch) //Also diagonal
				multi++;
			if (matrix[i][j] == gap2) //Also up
				multi++;
		}
    //Moving up
		if (i > 0 && path[i][j] == &matrix[i - 1][j]) {
			newSeq2.insert(0, 1, '-');
			newSeq1.insert(0, 1, seq1[i - 1]);
			i--;
      //Checking for additional paths
			if (matrix[i][j] == match || matrix[i][j] == mismatch) //Also diagonal
				multi++;
			if (matrix[i][j] == gap1)  //Also left
				multi++;
		}
		if (multi > 0)
			multiplePaths = true;
	}

	printIfMultiple(multiplePaths);

	alignFile.open("assignment1.o3");
	if (alignFile.is_open()) {
		alignFile << newSeq1 << "\n";
		alignFile << newSeq2 << "\n";
		alignFile.close();
	}
	else cout << "Unable to open file.\n";
}

//Indicates presence of multiple optimal alignments with YES or NO
void printIfMultiple(bool ifMulti) {
	ofstream multiFile;
	multiFile.open("assignment1.o4");
	if (multiFile.is_open()) {
		if (ifMulti)
			multiFile << "YES\n";
		else
			multiFile << "NO\n";
		multiFile.close();
	}
	else cout << "Unable to open file.\n";
}