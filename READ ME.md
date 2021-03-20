# READ ME - De Bruijn Assembler



## About the program

The aim of this program is to simulate a situation where different fragments of DNA must be reassembled without having a reference genome for comparison. The whole assembly must be conducted based on the overlap between the fragments.

The program is designed to split a DNA sequence into fragments of size chosen by the user. The fragments are reassembled using the De Bruijn Graph approach, that is based on the overlapping fragments.  An optional output containing a visual illustration of the De Bruijn Graph can be created.



## Installation

The program is written in Python, and uses the following versions:

- Python 3.8.5
- pip 21.0.1
- graphviz 0.16

### Install graphviz

The installation of *graphviz* is done using pip, so make sure you have pip installed before trying to install *graphviz*.

```
# check if pip is installed
pip help
# if the help option appears, continue to graphviz installation, otherwise install pip first.
pip install graphviz==0.16
```



## Running the program

The program can be run in two different ways. Using the graphical interface or through command line.

### Using the Graphical Interface

To start the Graphical Interface, you just need to run the program without any flag.

```
python bruijn_assembler.py
```

The contents of the initial window of the Graphical Interface are detailed below.

- **Sequence box:** text box where the sequence can be inserted.
- **Kmer box:** text box to insert the kmer size.
- **Checkbox** to create a De Bruijn Graph as png. If selected, after submitting, a save file window will appear, so you can choose the name and where to save your png image. The file opens automatically after created.
- **Open file button:** open a file selection window and retrieve the sequence from the selected file. The file must be in fasta format, otherwise it will not be read.
- **Submit button:** retrieve the information and continue for the result window.
- **Clear button:** clear both the sequence box and the kmer box.
- **Close button:** close the program.

### Using command line

To run the program via command line, you have to specify the flag **-c**.

```
python bruijn_assembler.py -c --file sequence.fasta --kmer 12 --draw_graph --out de_bruijn_graph.png
```

Flag options:

- **-c :** stops the program from opening the Graphical Interface.
- **--file :** retrieve the sequence from a file in fasta format.
- **--sequence :** input a sequence.
- **--kmer : **integer representing the size of the fragments.
- **--draw_graph :** creates a png file containing a De Bruijn Graph.
- **--out : **file name which the De Bruijn Graph will be saved.

When the program is used via command line, the flags **--file** or **--sequence**, and **--kmer** are mandatory. 

The flags **--draw_graph** and **--out** are optional, the default png file name is "result_de_bruijn_graph.png".



## Results

The result is divided in three parts: a list of all the possible kmers, both the original sequence and the result sequence (next to each other to comparison), and the similarity percentage between the original and result sequences. The result will be presented differently in the Graphical User Interface and the command line.

### Graphical interface

The list of kmers is found on the left side of the window, and both original and result sequences are in the main text box with a label differing them. Under the box with the sequences is possible to find the percentage of similarity between them.

A new button will be available, the New Sequence button that will go back to the initial window.

### Command line

Using the command line will generate the result in standard output, that can be easily saved to a file. The kmer list is printed first. The sequences are printed next to each other for comparison, after the kmer list. The last printed sentence contains the percentage of similarity between the sequences.



