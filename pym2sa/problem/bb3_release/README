******************************************************************************

               BAliBASE multiple alignment benchmark suite
                        (version 3.0, Dec 2004)

******************************************************************************


Please send bug reports, comments etc. to :-
        thompson@igbmc.u-strasbg.fr

******************************************************************************

The tests in the BAliBASE benchmark suite are divided into 5 different reference sets.

Reference 1 equi-distant sequences with 2 different levels of conservation,

Reference 2 families aligned with a highly divergent "orphan" sequence,

Reference 3 subgroups with <25% residue identity between groups,

Reference 4 sequences with N/C-terminal extensions,

Reference 5 internal insertions.

The different reference sets are organised in seperate directories:

RV11:	Reference 1, very divergent sequences (<20% identity)
RV12:	Reference 1, medium to divergent sequences (20-40% identity)
RV20:	Reference 2
RV30:	Reference 3
RV40:	Reference 4
RV50:	Reference 5

For each test, a number of files are provided. The files named BBnnnnn contain the 
full-length sequences, while the files named BBSnnnnn contain the sequences corresponding
to the homologous regions only. Different file formats are also provided:

BBnnnnn.xml	the XML file with the aligned full length sequences and the annotation
BBnnnnn.msf	the MSF file with the aligned full length sequences 
BBnnnnn.tfa	the TFA file with the non-aligned full length sequences 
BBSnnnnn.xml	the XML file with the aligned truncated sequences and the annotation
BBSnnnnn.msf	the MSF file with the aligned truncated sequences 
BBSnnnnn.tfa	the TFA file with the non-aligned truncated sequences 

NB. Reference 4 only contains the full-length sequences, because the test is designed to
evaluate the effect of N/C-terminal extensions!


A C program called bali_score is also provided that reads a BAliBASE reference file in
XML or MSF format and compares it to a user's alignment in MSF format. If the XML format 
file is used, the alignments are compared only in the core blocks defined in the XML file.
If the MSF reference file is used, the alignments are compared using all columns in the alignment.
Two scores are calculated: the sum-of-pairs and the total column score. (see reference 2 for
details).

References

1. Thompson JD, Plewniak F, Poch O. BAliBASE: a benchmark alignment database for the evaluation of multiple alignment programs.
Bioinformatics. 1999 Jan;15(1):87-8. 

2. Thompson JD, Plewniak F, Poch O. A comprehensive comparison of multiple sequence alignment programs.
Nucleic Acids Res. 1999 Jul 1;27(13):2682-90

3. Bahr A, Thompson JD, Thierry JC, Poch O. BAliBASE (Benchmark Alignment dataBASE): enhancements for repeats, transmembrane sequences and circular permutations.
Nucleic Acids Res. 2001 Jan 1;29(1):323-6. 
