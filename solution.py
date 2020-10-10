const MATCH


def buildMatrix(seq1, seq2):
  matrix = [[0]*len(seq2)] * seq1

def main():
  fileName = input("Enter File Name\n")
  file = open(fileName, 'r')
  nextIsSeq = False
  firstSeq = True
  seq1, seq2 = None, None
  for line in file:
    if nextIsSeq:
      if firstSeq: 
        seq1 = line.strip()
        firstSeq = False
      else: 
        seq2 = line.strip()
      nextIsSeq = False
    if line[0] == '>':
      nextIsSeq = True

  matrix = buildMatrix(seq1, seq2)

if __name__ == '__main__':
    main()