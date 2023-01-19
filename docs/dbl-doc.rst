
Overview
========

dmbiolib is a library of functions used in various bionformatics projects.

Source code:
 https://github.com/damienmarsic/dmbiolib

Python package:
 https://pypi.org/project/dmbiolib/


Installation
============

dmbiolib can be installed using pip::

    pip install dmbiolib

Note that dependencies might need to be installed individually.


Usage
=====

dmbiolib needs to be imported before its functions can be used.
Example::
    import dmbiolib as dbl
    print(dbl.transl('atgcgattcacg'))


Functions
=========

aln2seq(filename,type,full,reference)
*************************************
| filename: file containing multiple sequence alignment in caplib3 format
| type: dna or aa
| full: True if full sequences are to be returned (only valid if reference is provided)
| reference: name of file containing the reference sequence
| Converts alignments into caplib3 format (variable regions with identity displayed as .) into complete sequences.
| *Work in progress !*

check_file(filename,strict)
***************************
| filename: file name to be tested
| strict: True or False
| Checks whether a file (filename) exists in the working directory.
| Returns True if the file exists, False if the file does not exist and strict is False. Exits if the file does not exist and strict is True.

check_plot_format(x)
********************
| x: string to be tested
| Checks whether a string (x) is an accepted plot format.
| Returns x if x is svg, png, jpg, jpeg, pdf, ps, eps, pgf, raw, rgba, tif or tiff.
| Returns empty string if x is 'Single multipage pdf'.
| Exits in all other cases.

check_read_file(filename)
*************************
| filename: name of file to be tested for containing sequencing reads
| Returns a fail message if filename if found not to be a valid read file. Returns an empty string otherwise.

check_seq(sequence,type,required)
*********************************
| sequence: amino acid or nucleotide sequence
| type: dna ('atgc'), ambiguous ('ryswkmbdhvn'), aa ('ARNDCQEGHILKMFPSTWYV'), or any string (case insensitive)
| required: string of characters (or type name as above), at least one of which must exist in the sequence
| Checks whether the sequence is of the correct type and if at least one of the required characters is present in the sequence.
| Returns a pair of booleans: x,y
| x: True if all characters in sequence are in type, False if any character is not in type
| y: True if at least 1 character in sequence is in required, False otherwise

Examples::
    import dmbiolib as dbl
    print(dbl.check_seq('cgttcgaac',dbl.dna,dbl.dna))
    True, True
    print(dbl.check_seq('cgttnnaac',dbl.dna,dbl.dna))
    False, True
    print(dbl.check_seq('cgttnnaac',dbl.dna,dbl.ambiguous))
    True, True


check_sync(read1,read2)
***********************
| read1, read2: nucleotide sequences
| Checks whether the 2 read files (Illumina paired-ends) are synchronized (reads in the same file location belongs to the same pair).
| Returns a fail message if the files are not synchronized. Returns an empty string otherwise.

complexity(sequence)
********************
| sequence: nucleotide sequence (including ambiguous nucleotides) to be translated (in frame)
| Returns a list of dictionaries. Each list item corresponds to a nucleotide triplet from the sequence. Each dictionary lists amino acids corresponding to the triplet translation, with the number of different codons for each amino acid.

Example::
   import dmbiolib as dbl
   x=dbl.complexity('atgdbctss')
   for n in x:
       print(n)
   defaultdict(<class 'int'>, {'M': 1})
   defaultdict(<class 'int'>, {'F': 1, 'C': 1, 'S': 2, 'V': 1, 'G': 1, 'A': 1, 'I': 1, 'T': 1})
   defaultdict(<class 'int'>, {'W': 1, 'C': 1, 'S': 2})


compress(sequence):
*******************
| sequence: nucleotide sequence
| Returns a "compressed" sequence in which all homopolymers (but only if a, g, c or t) are shortened to just one copy.

Example::
   import dmbiolib as dbl
   print(dbl.compress('gggcaatccccnnnncaagtt'))
   gcatcnnnncagt
   
conf_start(filename,title)
**************************
| Creates a configuration file, using filename (name of configuration file) and title (text to be included in the title at the beginning of the file).
| Returns the file handle, the current directory name and a list of detected read files or read file pairs preceded by a file prefix.

conf_end(file_handle,filename,title)
************************************
| Completes writing the configuration file.

csv_read(filename,dic,header)
*****************************
| filename: name of csv file to be read
| dic (True/False): whether to store the contents of the csv file in a dictionary (True) or a lst (False).
| header (True/False): whether the file starts with a header or not (or directly with the data)
| Opens a csv file and stores its content into a dictionary, while converting numbers to integers or floats as appropriate.

csv_write(filename,keys,list_or_dic,header,description,file_handle)
*******************************************************************
| filename: name of csv file to be created
| keys: optional first column (if not already part of the list or dictionary)
| list_or_dic: list (or tuple) or dictionary containing the data to be written into the csv file
| header: optional top row to be written before the main data
| description: file description to be used in the message confirming completion of csv file
| file_handle: file_handle of the report file (or None if no report file)
| Creates a csv file from the arguments.

diff(sequences)
***************
| sequences: list of sequences
| Returns the smallest number of differences between any 2 sequences from the list. This is useful to evaluate a list of barcodes for example, to make sure all barcodes differ from each other by at least some number of differences. Note that all sequences must be of the same length.

Examples::
   import dmbiolib as dbl
   print(dbl.diff(['agct','gatc','ctga','tcag']))
   4
   print(dbl.diff(['agct','gatc','ctga','aata']))
   2

dirname()
*******
| Returns the name (not the full path) of the current directory.

Example, if current directory is /home/someuser/somedir::
   print(dirname())
   somedir

entropy(matrix)
***************
| matrix: list of lists of values
| Returns the Shannon entropy of the matrix.

exprange(a,b,c)
***************
| a,b: range boundaries
| c: multiplying factor
| Returns an exponential range as a generator.

Example::
   import dmbiolib as dbl
   x=dbl.exprange(1,100,3)
   for n in x:
       print(n)
   1
   3
   9
   27
   81

find_ambiguous(seq)
*******************
| seq: nucleotide sequence (containing ambiguous nucleotides)
| Identifies location of all ambiguous stretches and their length, which it returns as a dictionay.

Example::
   import dmbiolib as dbl
   seq='gatcgatcgtnnnnngactgavvmttcgsbynccgtcga'
   print(dbl.find_ambiguous(seq))
   {10: 5, 21: 3, 28: 4}

find_read_files()
*****************
| Looks for read files (gzipped only) in the current directory.
| Returns a list in wich each item is a string containing a prefix followed by either a single read file or a pair (in case of paired ends sequencing), separated by a space.

format_dna(seq,margin,cpl,cpn)
******************************
| seq: raw nucleotide sequence
| margin: left margin
| cpl: number of characters per line
| cpn: number of characters per number
| Returns formatted nucleotide sequence.

Example::
   seq='gatcgatcgatcgatcgtacgtatcgatcgatcgatcgatcgactgatcagctacgatcgatcgatcgatgtgacccccttagc'
   print(dbl.format_dna(seq,5,30,10))
                10        20        30
        gatcgatcgatcgatcgtacgtatcgatcg
                40        50        60
        atcgatcgatcgactgatcagctacgatcg
                70        80
        atcgatcgatgtgacccccttagc

fsize(filename)
***************
| Returns the size in bytes of the file named filename.

getfasta(fname,type,required,multi)
***********************************
| fname: name of the fasta file to be opened
| type: dna or aa
| required: same as type, or 'ambiguous' if some ambiguous nucleotides must be present
| multi: Whether the file contains multiple sequences (True) or a single one (False).
| Returns a dictionary of all sequences identified (keys: sequence names, values: sequences) and a string containing possible fail messages.

getread(f,y,counter)
********************
| f: file handle
| y: number of lines per sequence (or 0 if variable number)
| counter: number of reads already processed
| Reads next read and determine read name and sequence.
| Returns read sequence, file handle, updated counter, read name.

initreadfile(rfile)
*******************
| rfile: read file (can be fasta or fastq, uncompressed or gzipped)
| Opens and checks the file. Detects if the format is fastq (new sequence every 4 lines), single line fasta (new sequence every 2 lines) or multiline fasta (new sequence every unknown number of lines).
| Returns file handle and number of lines for each sequence (or 0 if format is multiline fasta).

intorfloat(x)
*************
| x: string to be tested whether it can be converted into an integer or a float
| Returns 'int' if x can be converted to an integer, 'float' if can be converted into a float, 'other' in all other cases.

lncount(f)
**********
| f: file handle
| Returns the number of lines in the file (works fast with large files).

match(seq1, seq2)
*****************
| seq1, seq2: nucleotide sequences (with or without ambiguous nucleotides)
| Checks if the 2 sequences match at each position (see nt_match() below).
| Returns True if the sequences match, False otherwise (or if sequence lengths are different).

Examples::
   import dmbiolib as dbl
   dbl.match('acgatcg','accatcg')
   False
   dbl.match('acgatcg','acsancg')
   True

mean(x)
*******
| x: list or tuple of numerical values
| Returns the mean (sum of all values divided by number of values).

Example::
   import dmbiolib as dbl
   print(dbl.mean([12,30,24]))
   22.0

nt_match(nt1, nt2)
******************
| nt1, nt2: nucleotide (a, g, c, t or ambiguous)
| Returns True if the 2 nucleotides match, False otherwise.
| Matching means identity for a, t, g and c, and compatibility for ambiguous nucleotides.

Examples::
   import dmbiolib as dbl
   dbl.nt_match('a','a')
   True
   dbl.nt_match('a','g')
   False
   dbl.nt_match('n','a')
   True
   dbl.nt_match('s','n')
   True
   dbl.nt_match('r','y')
   False
   dbl.nt_match('g','s')
   True

open_read_file(filename)
************************
| filename: name of the read file to be opened
| Opens a read file (either uncompressed or gzipped) and returns the file handle.

plot_end(fig,name,format,mppdf)
*******************************
| fig: figure handle
| name: file name without extension (if each figure is saved individually)
| format: extension corresponding to the chosen figure format (if each figure is saved individually)
| mppdf: PdfPages handle (if all figures saved in single file pdf)
| Completes the plotting process.

plot_start(x,y,z)
*****************
| x: color map to be used
| y: number of colors needed
| z: plot title
| Initializes the plot
| Returns list of colors and figure handle

pr2(f,text)
***********
| f: file handle
| text: text to be printed
| Prints a text simultaneously to the screen and to a file (adds '\n' when printing to file).

prefix(x)
*********
| x: list of file names
| Returns a list of numbers, with each number being the suggested slice (from left end) of the corresponding file name to be used as a prefix.

Example::
   import dmbiolib as dbl
   x=['P0-left_L4_2.fq.gz', 'P0-right_L4_2.fq.gz', 'P1-left_L4_2.fq.gz', 'P1-right_L4_2.fq.gz', 'P2-left_L4_2.fq.gz', 'P2-right_L4_2.fq.gz']
   print(dbl.prefix(x))
   [7, 8, 7, 8, 7, 8]

progress_check(c,show,text)
***************************
| c: read counter
| show: dictionary of read numbers that trigger a new % value to the progress counter
| text: text describing the process (should be the same as in progress_start(nr,text))
| Updates the progress counter that was created by progress_start(nr,text).

progress_end()
**************
| Prints the final 100.0% when the process has been completed.

progress_start(nr,text)
***********************
| nr: number of reads
| text: text describing the process
| Starts a progress counter (from 0.0% to 100.0%) of going through a read file.
| Returns a dictionary of read numbers and % completion (only the read numbers that will trigger an update to the counter).

readcount(R,fail)
*****************
| R: name of read file
| fail: fail message
| Counts number of reads in a read file (can be fasta or fastq format, either uncompressed of gzipped). Add a fail text to the fail variable if the file if detected as not being a read file.
| Returns number of reads and updated fail message.

rename(filename)
****************
| filename: name of the file to be renamed
| If the file exists and has non zero size, it is renamed by appending a unique number to it, so a new file with the name filename can be created.

revcomp(seq)
************
| seq: nucleotide sequence
| Returns the reverse-complement.

Example::
   revcomp('agctgctaa')
   ttagcagct

shortest_probe(seqs,lim,host,t)
*******************************
| seqs: list of nucleotide sequences
| lim: minimum probe size
| host: host genome
| t: description
| Returns shortest probe size allowing to identify all sequences and with probe sequence not present in the host genome.

sortfiles(l,str)
****************
| l: list of file names to be sorted
| str: string before which file names will be sorted
| Returns a list of sorted file names. Sorting is based on numbers if numbers are present in the file names.

transl(seq)
***********
| seq: nucleotide sequence
| Returns amino acid sequence translation of the nucleotide sequence.

Example::
   transl('atgctgaaagcc')
   MLKA


