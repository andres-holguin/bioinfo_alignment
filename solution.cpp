#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <iomanip>
using namespace std;

/*************SCORING FUNCTION:**********************/
/*2 for match, -1 for mismatch, -2 match against gap*/
/****************************************************/
struct Scores {
	static int const 
		MATCH = 2,
		MISMATCH = -1,
		GAP = -2;
};
enum Direction {	UP, DIAG, LEFT, NONE };
struct AlignmentScore {int val; Direction source = NONE;};
class AlignmentScoreMatrix {
	string SEQ1, SEQ2;
	u_int M, N;
	vector< vector<AlignmentScore> > matrix;
public:
	AlignmentScoreMatrix(string const &SEQ1, string const &SEQ2) {
		this->SEQ1 = SEQ1, this->SEQ2 = SEQ2;
		M = SEQ1.size()+1, N = SEQ2.size()+1;
		matrix = vector< vector<AlignmentScore> >(M, vector<AlignmentScore>(N)); 
	}

	AlignmentScore* at(u_int i, u_int j) {
		return &matrix[i][j];
	}

	u_int const getNumRows() {
		return M;
	}
	u_int const getNumCols() {
		return N;
	}
	string const getSEQ1() {
		return SEQ1;
	}
	string const getSEQ2() {
		return SEQ2;
	}
};

vector<string> extractSeq(string const &fileName);
AlignmentScoreMatrix getScores(string const &SEQ1, string const &SEQ2);
void solveFromSubProblem(AlignmentScoreMatrix &matrix, u_int i, u_int j, Direction dir, int &maxScore);

void writeScoreMatrix(AlignmentScoreMatrix &matrix, string const &fileName);

// vector< vector<int> > scoreMatrix(string seq1, string seq2);
// void printIfMultiple(bool ifMulti);
// void printOptAlignment(vector<vector<int>>& matrix, vector<vector<int*>>& path, string seq1, string seq2);
// EC: Write out all possible optimal alignments


int main(int argc, char** argv) {
	/* Determine file to open */
	string fileName;
	if (argc > 1) { // File name is provided as argument
		fileName = argv[1];
	} else { // Prompt for the file's name
		cout << "Enter File Name" << endl;
		getline(cin, fileName);
	}

	/* Open the file and extract sequences */
	cout << "Extracting sequences from " << fileName << endl;
	vector<string> sequences = extractSeq(fileName);
	if ( !sequences.size() ) return -1; // No Sequences were extracted;
	// Assign sequences
	const string SEQ1 = sequences[0], SEQ2 = sequences[1];

	/* Get score matrix */
	cout << "\nComparing two sequences:" << endl;
	cout << SEQ1 << endl;
	cout << "and" << endl;
	cout << SEQ2 << endl;
	cout << "With scores: " << endl;
	cout << "\tMatch: " << Scores::MATCH << endl;
	cout << "\tMismatch: " << Scores::MISMATCH << endl;
	cout << "\tGap: " << Scores::GAP << endl;
\
	AlignmentScoreMatrix matrix = getScores(SEQ1, SEQ2);

	/* Output results */
	// (a) Write the optimal score

	// (b) Write the dynamic programming matrix
	writeScoreMatrix(matrix, "assignment1.o2");

}

vector<string> extractSeq(const string &fileName) {
	ifstream inFile;
	string inLine;
	vector<string> sequences;

	inFile.open(fileName);
	if (!inFile.is_open()) cout << "Unable to open file: " << fileName << endl;

	while( getline(inFile, inLine) ) {
		if ( inLine.size() && inLine.at(0) == '>' ) {
			getline(inFile, inLine);
			sequences.push_back(inLine);
		}
	}
	return sequences;
}

/** 
 * Builds the solution matrix of the form:
 * x _ S E Q 2
 * _ 0 0 0 0 0
 * S 0 0 0 0 0
 * E 0 0 0 0 0
 * Q 0 0 0 0 0
 * 1 0 0 0 0 0
 */
AlignmentScoreMatrix getScores(string const &SEQ1, string const &SEQ2) {
	// Init constants to sequence sizes
	u_int const M = SEQ1.size()+1, N = SEQ2.size()+1;
	// Create matrix of size ROW x COL = M X N
	AlignmentScoreMatrix matrix(SEQ1, SEQ2);

	auto isWithinMatrix = [&](u_int const &i, u_int const &j) -> bool { 
		return (i <= M && j <= N); // Since unsigned, only need to check if greater than.
	};
	
	for (u_int i = 0; i < M; i++) {
		for (u_int j = 0; j < N; j++) {
			int maxScore = INT_MIN;
			u_int seq1_index, seq2_index;

			// Recurse from top-left (DIAG) for MATCH or MISMATCH
			solveFromSubProblem(matrix, i, j, DIAG, maxScore);

			// Recurse from above (UP) for GAP on SEQ2
			solveFromSubProblem(matrix, i, j, UP, maxScore);

			// Recurse from LEFT for GAP or SEQ1
			solveFromSubProblem(matrix, i, j, LEFT, maxScore);
		}
	}

	return matrix;
}

void solveFromSubProblem(AlignmentScoreMatrix &matrix, u_int i, u_int j, Direction dir, int &maxScore) {
	u_int const M = matrix.getNumRows(), N = matrix.getNumCols();
	int newScore, scoreMod;
	u_int seq1_index, seq2_index;
	// Lambda to determine if indices within matrix
	auto isWithinMatrix = [&](u_int const &i, u_int const &j) -> bool { 
		return (i <= M && j <= N); // Since unsigned, only need to check if greater than.
	};

	switch (dir) {
	case DIAG: // Recurse from top-left for MATCH or MISMATCH
		seq1_index = i-1, seq2_index = j-1;
		scoreMod = isWithinMatrix(seq1_index, seq2_index) && matrix.getSEQ1()[seq1_index] == matrix.getSEQ2()[seq2_index] ? 
								Scores::MATCH :
								Scores::MISMATCH ;
		break;
	case UP: // Recurse from above for GAP on SEQ2
		seq1_index = i, seq2_index = j-1;	
		scoreMod = Scores::GAP;
		break;
	case LEFT: // Recurse from LEFT for GAP on SEQ1
		seq1_index = i-1, seq2_index = j;
		scoreMod = Scores::GAP;		
		break;
	default:
		return;
	}

	// Set as new max
	if ( isWithinMatrix(seq1_index, seq2_index) ) {
		newScore = matrix.at(seq1_index, seq2_index)->val + scoreMod;
		if (newScore > maxScore) {
			matrix.at(i,j)->val = maxScore = newScore;
			matrix.at(i,j)->source = dir;
		}
	}	
}

char interpretDirection(Direction dir) {
	switch (dir) {
		case UP:
			return '^';
		case DIAG:
			return '\\';
		case LEFT:
			return '<';
		default:
			return '.';
	}
}

void writeScoreMatrix(AlignmentScoreMatrix &matrix, string const &fileName) {
	ofstream outfile(fileName);
	if (outfile.is_open())
	{
		cout << "\nScore Matrix:" << endl;
		for (u_int i=0; i < matrix.getNumRows(); i++) {
			for (u_int j=0; j < matrix.getNumCols(); j++) {
				outfile << matrix.at(i,j)->val << ' ';
				cout << right << setw(4) << matrix.at(i,j)->val << interpretDirection(matrix.at(i,j)->source) << ' ';
			}
			outfile << endl;
			cout << endl;
		}
		outfile.close();
	}
	else cout << "Unable to output the score matrix to " << fileName << endl;
}

/*
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
*/
