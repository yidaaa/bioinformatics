INTRODUCTION:

Sequence alignment handles 2 sequences and returns their best alignment between prefixes.
In other words, ending spaces for both sequences are forgiven.

Three programs are implemented in python: 
  - standard dynamic programming (align_DP.py)
  - the banded dynamic programming (align_band)
  - X-drop algorithm (align_drop)

All 3 programs can:
1. Align sequences over any alphabet. The alphabet is specified in the parameter file.
2. Use any score matrix. The score matrix is specified in the parameter file.
3. Handle affline gap penalty. (For linear gap penalty, just set opening gaps to 0)

Each program, outputs the alignment score, the alignment and the number of entries filled in during matrix calculation.
Output will be printed into the output file as well as your console.


INSTRUCTIONS TO RUN THE PROGRAM:

Three text files should be present in the same directory as the program files to successfully run the program:
- parameter file
- input file
- output file

The files should follow the below formats:

1.  parameterfile.txt
Do not remove any line or semi-colons. Replace the numbers to change the parameters in the example below.

        EXAMPLE
        ------------------------------------------------------
        1; the threshold X for X-drop
        3; bandwidth B
        -1; score for initiating a gap
        -1; score for each base insert/delete

        ; Below are alphabets in DNA
        a c g t

        ; Below is the similarity matrix for the alphabets
        2 -1 -1 -1
        -1 2 -1 -1
        -1 -1 2 -1
        -1 -1 -1 2
        ------------------------------------------------------

2.  inputfile.txt
    The input file should be in standard FASTA format. 2 sequences' names are indicated by ">"
    Respective sequences to be aligned are below the name:

        EXAMPLE
        ------------------------------------------------------
        >seq1
        gtagcatggt

        >seq2
        gtacaatccc
        ------------------------------------------------------

3.  outputfile.txt
    An empty text document. Output will be printed here after running the program.
     
Example files of the above parameter and input text files are avaliable.

After running the program, files should be input in the following format:
> (parameter file name).txt (input file name).txt (output file name).txt

        EXAMPLE ON COMMAND LINE:
        ------------------------------------------------------
        >align_DP.py
        >parameter1.txt input1.txt output1.txt
        ------------------------------------------------------
        
Your output text file will show the results of the alignment:

        EXAMPLE
        ------------------------------------------------------
        score =  8
        entries =  341

        >seq1
        gtagca

        >seq2
        gta-ca
        ------------------------------------------------------
