
Overview
========

dmbiolib is a library of functions used in various bionformatics projects.

Source code:
 https://github.com/damienmarsic/dmbiolib

Python package:
 https://pypi.org/project/dmbiolib/

Bioconda package:
 https://bioconda.github.io/recipes/dmbiolib/README.html

Bug report / feature requests:
 https://github.com/damienmarsic/dmbiolib/issues/new/choose


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


Latest news
===========

| 2023-05-16: update to 0.4.0. Improved conf_start, conf_end, cvs_write, prefix. Added aa_dist, detect_vr, findall, frame, mut_per_read, prod, seq_clust_card_dist, seq_write, size_dist.
| 2023-05-16: update to 0.4.1. Minor fixes.


Functions
=========

aa_dist(seqs,parvrs,fname,r)
****************************
* seqs: dictionary of dictionaries {vr1:{seq1:n1, seq2:n2, ...}, vr2:{...}, ...} where vr1, vr2 are variable region names, seq1, seq2 are amino acid sequences, n1, n2 are read counts
* parvrs: dictionary of tuples {vr1:(seq1, pos1), vr2:(seq2, pos2), ...} where vr1, vr2 are variable region names (must be the same as in seqs), seq1, seq2 are parental amino acid sequences, pos1, pos2 are position numbers in the parental protein chain (can be None of no parental sequence)
* fname: name of output file (can be None if no need to save results)
* r: handle of report file (can be None if not used)

| Creates an amino acid distribution dictionary and csv file, including each position in each variable region.

| Also creates a mutation distribution dictionary and csv file if parvrs is provided.

| Returns amino acid distribution and mutation distribution dictionaries.


aln2seq(filename,type,full,reference)
*************************************
* filename: file containing multiple sequence alignment in caplib3 format
* type: dna or aa
* full: True if full sequences are to be returned (only valid if reference is provided)
* reference: name of file containing the reference sequence

| Converts alignments into caplib3 format (variable regions with identity displayed as .) into complete sequences.

| *Work in progress !*

check_file(filename,strict)
***************************
* filename: file name to be tested
* strict: True or False

| Checks whether a file (filename) exists in the working directory.

| Returns True if the file exists, False if the file does not exist and strict is False. Exits if the file does not exist and strict is True.

check_plot_format(x)
********************
* x: string to be tested

| Checks whether a string (x) is an accepted plot format.

| Returns x if x is svg, png, jpg, jpeg, pdf, ps, eps, pgf, raw, rgba, tif or tiff. Returns empty string if x is 'Single multipage pdf'. Exits in all other cases.

check_read_file(filename)
*************************
* filename: name of file to be tested for containing sequencing reads

| Returns a fail message if filename if found not to be a valid read file. Returns an empty string otherwise.

check_seq(sequence,type,required)
*********************************
* sequence: amino acid or nucleotide sequence
* type: dna ('atgc'), ambiguous ('ryswkmbdhvn'), aa ('ARNDCQEGHILKMFPSTWYV'), or any string (case insensitive)
* required: string of characters (or type name as above), at least one of which must exist in the sequence

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
* read1, read2: nucleotide sequences

| Checks whether the 2 read files (Illumina paired-ends) are synchronized (reads in the same file location belongs to the same pair).

| Returns a fail message if the files are not synchronized. Returns an empty string otherwise.

complexity(sequence)
********************
* sequence: nucleotide sequence (including ambiguous nucleotides) to be translated (in frame)

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
* sequence: nucleotide sequence

| Returns a "compressed" sequence in which all homopolymers (but only if a, g, c or t) are shortened to just one copy.
Example::

   import dmbiolib as dbl
   print(dbl.compress('gggcaatccccnnnncaagtt'))
   gcatcnnnncagt
   
conf_start(title)
**************************
| Creates a configuration file, using title (text to be included in the title at the beginning of the file).

| Returns a string with the preliminary content of the future configuration file, the current directory name and a list of detected read files or read file pairs preceded by a file prefix.

conf_end(filename,content,title)
************************************
| Completes writing content into the configuration file.

csv_read(filename,dic,header)
*****************************
* filename: name of csv file to be read
* dic (True/False): whether to store the contents of the csv file in a dictionary (True) or a lst (False).
* header (True/False): whether the file starts with a header or not (or directly with the data)

| Opens a csv file and stores its content into a dictionary, while converting numbers to integers or floats as appropriate.

csv_write(filename,keys,list_or_dic,header,description,file_handle)
*******************************************************************
* filename: name of csv file to be created
* keys: optional first column (if not already part of the list or dictionary)
* list_or_dic: list (or tuple) or dictionary containing the data (which can be strings, lists, tuples or dictionaries) to be written into the csv file
* header: optional top row to be written before the main data
* description: file description to be used in the message confirming completion of csv file
* file_handle: file_handle of the report file (or None if no report file)

| Creates a csv file from the arguments.

detect_vr(libnt,mindist)
************************
* libnt: in frame, protein-coding library nucleotide sequence (containing ambiguous positions)
* mindist: minimum distance between 2 variable regions

| Detects variable regions (as strings of codons) from a library nucleotide sequence.

| Returns a dictionary of lists: {vr1:[left_seq,vr_seq,right_seq], vr2:...} where vr1, vr2: variable region names, left_seq: nucleotide sequence upstream the variable region, vr_seq, variable region sequence, right_seq: nucleotide sequence downstream the variable region.

diff(sequences)
***************
* sequences: list of sequences

| Returns the smallest number of differences between any 2 sequences from the list. This is useful to evaluate a list of barcodes for example, to make sure all barcodes differ from each other by at least some number of differences. Note that all sequences must be of the same length.
Examples::

   import dmbiolib as dbl
   print(dbl.diff(['agct','gatc','ctga','tcag']))
   4
   print(dbl.diff(['agct','gatc','ctga','aata']))
   2

dirname()
*********
| Returns the name (not the full path) of the current directory.
Example, if current directory is /home/someuser/somedir::

   print(dirname())
   somedir

entropy(matrix)
***************
* matrix: list of lists of values

| Returns the Shannon entropy of the matrix.

exprange(a,b,c)
***************
* a,b: range boundaries
* c: multiplying factor

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
* seq: nucleotide sequence (containing ambiguous nucleotides)

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

findall(probe,seq,start,end,overlap=False)
******************************************
* probe: string, occurrences of which are searched in seq
* seq: string in which probe is searched
* start: seq start index of search (0 if no limit)
* end: seq end index of search (None of no limit)
* overlap: optional argument to allow overlaps (default: False)

| Finds all occurrences of a string (probe) in a bigger string (seq), between seq start and end, with overlaps included optionally.

| Returns an iterator (must be converted to list if a list is needed).

format_dna(seq,margin,cpl,cpn)
******************************
* seq: raw nucleotide sequence
* margin: left margin
* cpl: number of characters per line
* cpn: number of characters per number

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

frame(seq,strict=False)
***********************
* seq: nucleotide sequences to be examined
* strict: when True, will return None if the guess is too speculative (optional argument, default: False)

| Guesses the reading frame of seq.

| Returns the frame as 0, 1 or 2, or None if could not be guessed.

fsize(filename)
***************
| Returns the size in bytes of the file named filename.

getfasta(fname,type,required,multi)
***********************************
* fname: name of the fasta file to be opened
* type: dna or aa
* required: same as type, or 'ambiguous' if some ambiguous nucleotides must be present
* multi: Whether the file contains multiple sequences (True) or a single one (False).

| Returns a dictionary of all sequences identified (keys: sequence names, values: sequences) and a string containing possible fail messages.

getread(f,y,counter)
********************
* f: file handle
* y: number of lines per sequence (or 0 if variable number)
* counter: number of reads already processed

| Reads next read and determine read name and sequence.

| Returns read sequence, file handle, updated counter, read name.

initreadfile(rfile)
*******************
* rfile: read file (can be fasta or fastq, uncompressed or gzipped)

| Opens and checks the file. Detects if the format is fastq (new sequence every 4 lines), single line fasta (new sequence every 2 lines) or multiline fasta (new sequence every unknown number of lines).

| Returns file handle and number of lines for each sequence (or 0 if format is multiline fasta).

intorfloat(x)
*************
* x: string to be tested whether it can be converted into an integer or a float

| Returns 'int' if x can be converted to an integer, 'float' if can be converted into a float, 'other' in all other cases.

match(seq1, seq2)
*****************
* seq1, seq2: nucleotide sequences (with or without ambiguous nucleotides)

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
* x: list or tuple of numerical values

| Returns the mean (sum of all values divided by number of values).
Example::

   import dmbiolib as dbl
   print(dbl.mean([12,30,24]))
   22.0

mut_per_read(seqs,parseq,fname,r)
*********************************
* seqs: dictionary {seq1:n1, seq2:n2, ...} where seq1, seq2: amino acid sequences, n1, n2: numbers of reads
* parseq: parental sequence (must be same length as sequences in seqs)
* fname: name of output file (can be None if no need to save results)
* r: handle of report file (can be None if not used)

| Creates a dictionary and csv file of distribution of number of mutations per read.

| Returns a dictionary {n1:m1, n2:m2, ...} where n1, n2: number of mutations, m1, m2: number of reads.

nt_match(nt1, nt2)
******************
* nt1, nt2: nucleotide (a, g, c, t or ambiguous)

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

plot_end(fig,name,format,mppdf)
*******************************
* fig: figure handle
* name: file name without extension (if each figure is saved individually)
* format: extension corresponding to the chosen figure format (if each figure is saved individually)
* mppdf: PdfPages handle (if all figures saved in single file pdf)

| Completes the plotting process.

plot_start(x,y,z)
*****************
* x: color map to be used
* y: number of colors needed
* z: plot title

| Initializes the plot

| Returns list of colors and figure handle

pr2(f,text)
***********
* f: file handle
* text: text to be printed

| Prints a text simultaneously to the screen and to a file (adds '\n' when printing to file).

prefix(x)
*********
* x: list of file names

| Returns a list of unique prefixes corresponding to the file names.
Example::

   import dmbiolib as dbl
   x=['P0-left_L4_2.fq.gz', 'P0-right_L4_2.fq.gz', 'P1-left_L4_2.fq.gz', 'P1-right_L4_2.fq.gz', 'P2-left_L4_2.fq.gz', 'P2-right_L4_2.fq.gz']
   print(dbl.prefix(x))
   ['P0-left', 'P0-right', 'P1-left', 'P1-right', 'P2-left', 'P2-right']

prod(x)
*******
* x: list or tuple of numbers

| Returns the product of all numbers in x

progress_check(c,show,text)
***************************
* c: read counter
* show: dictionary of read numbers that trigger a new % value to the progress counter
* text: text describing the process (should be the same as in progress_start(nr,text))

| Updates the progress counter that was created by progress_start(nr,text).

progress_end()
**************
| Prints the final 100.0% when the process has been completed.

progress_start(nr,text)
***********************
* nr: number of reads
* text: text describing the process

| Starts a progress counter (from 0.0% to 100.0%) of going through a read file.

| Returns a dictionary of read numbers and % completion (only the read numbers that will trigger an update to the counter).

readcount(R)
*****************
* R: name of read file
* fail: fail message

| Counts number of reads in a read file (can be fasta or fastq format, either uncompressed of gzipped). Add a fail text to the fail variable if the file if detected as not being a read file.

| Returns number of reads and updated fail message.

rename(filename)
****************
* filename: name of the file to be renamed

| If the file exists and has non zero size, it is renamed by appending a unique number to it, so a new file with the name filename can be created.

revcomp(seq)
************
* seq: nucleotide sequence

| Returns the reverse-complement.
Example::

   revcomp('agctgctaa')
   ttagcagct

rfile_create(filename)
************************
* filename: name of the read file to be created

| Creates a read file (either uncompressed or gzipped if .gz suffix is used) and returns the file handle.

rfile_open(filename)
********************
* filename: name of the read file to be opened

| Opens a read file (either uncompressed or gzipped) and returns the file handle.

seq_clust_card_dist(seqs,fname,r)
*********************************
* seqs: either a list [n1, n2, ...] or a dictionary {seq1:n1, seq2:n2, ...} where seq1, seq2: amino acid sequences, n1, n2: numbers of reads
* fname: name of output file (can be None if no need to save results)
* r: handle of report file (can be None if not used)

| Creates a dictionary and csv file of sequence cluster cardinality distribution.

| Returns a dictionary {n1:f1, n2:f2, ...} where n1, n2: cardinality, f1, f2: fraction of sequences.

seq_write(fname,top,seqs,dic,descr,r)
*************************************
* fname: name of file to be created
* top: string to be added to top of file
* seqs: list of sequences (or None)
* dic: dictionary of sequences with their read numbers {seq1:n1, seq2:n2, ...} (or None)
* descr: description to be included in message informing of task completion
* r: handle of report file (can be None if not used)

| Writes the sequences to a file, with an optional header string. Sequences can be either from a list or from a dictionary and are followed by read numbers if present in the dictionary.

shortest_probe(seqs,lim,host,t)
*******************************
* seqs: list of nucleotide sequences
* lim: minimum probe size
* host: host genome
* t: description

| Returns shortest probe size allowing to identify all sequences and with probe sequence not present in the host genome.

size_dist(seqs,fname,r)
***********************
* seqs: dictionary {seq1:n1, seq2:n2, ...} where seq1, seq2: amino acid sequences, n1, n2: numbers of reads
* fname: name of output file (can be None if no need to save results)
* r: handle of report file (can be None if not used)

| Creates a distribution of sequence lengths as a dictionary and a csv file.

| Returns a list of lengths (sorted by read numbers) and a dictionary of lengths and numbers of reads {l1:n1, l2:n2, ...} (l1, l2: sequence lengths, n1, n2: numbers of reads).

sortfiles(l,str)
****************
* l: list of file names to be sorted
* str: string before which file names will be sorted

| Returns a list of sorted file names. Sorting is based on numbers if numbers are present in the file names.

transl(seq)
***********
* seq: nucleotide sequence

| Returns amino acid sequence translation of the nucleotide sequence.
Example::

   transl('atgctgaaagcc')
   MLKA

xcount(f,x)
***********
* f: file handle (file must be opened in binary mode)
* x: string to be counted

| Returns the number of instances of x in the file (useful to count lines or reads in large files).



