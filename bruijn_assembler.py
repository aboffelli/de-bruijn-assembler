#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: De Bruijn Assembler

Created on: Mon Mar 08, 2021
Author: Arthur Boffelli Castro
"""
from graphviz import Digraph
import os
import argparse
import tkinter as tk
from tkinter import filedialog
import random

################################################################################
# Argparse set up.

usage = '''This program is designed to split a DNA sequence into fragments of 
size chosen by the user (kmer). The fragments are reassembled using the De 
Bruijn Graph approach, that is based on the overlapping fragments.
To run the program using the Graphical User Interface use: python 
bruijn_assembler.py.
To run the program using the command line use -c. The arguments --file or 
--sequence, and --kmer are mandatory.
The arguments --draw_graph and --out are optional.'''
'''
parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-c', dest='command',
                    action='store_true',
                    help='Run the program without opening the user interface'
                    )

parser.add_argument(
    '--sequence', type=str, metavar='', dest='input_sequence',
    help='Input sequence to be split'
)

parser.add_argument(
    '--file', type=str, dest='infile',
    help='Input file containing a sequence line in fasta format'
)

parser.add_argument(
    '--kmer', type=int, metavar='INT', dest='kmer',
    help='Length of the kmers to split the sequence (Integer)'
)

parser.add_argument(
    '--draw_graph', dest='draw_graph', action='store_true',
    help='Draw a De Bruijn Graph based on the kmer selected'
)

parser.add_argument(
    '--out', dest='outfile', metavar='OUTFILE', type=str,
    help='Path to save the pdf image if --draw_graph was selected'
)
args = parser.parse_args()


################################################################################
# Assemble functions


def kmer_breaks(sequence, k):
    """
    Function that breaks a sequence into a selected number of kmers.

    :param sequence: DNA sequence (str)
    :param k: length of the kmers (int)

    :return: sorted dictionary containing all the possible kmers in the sequence
    and the frequency that each kmer appears.
    """
    kmers_freq = {}
    sequence = sequence.upper()
    # loop through the sequence (-k + 1 to avoid index error).
    for i in range(len(sequence) - k + 1):

        fragment = sequence[i:i + k]
        if fragment not in kmers_freq:
            kmers_freq[fragment] = 1
        else:
            kmers_freq[fragment] += 1

    # Change the sequence of the kmers
    kmers_freq = dict(sorted(kmers_freq.items()))
    return kmers_freq


def reassembler(sequence, k):
    """
    Function to reassemble a split sequence based on the De Bruijn graph theory.

    :param sequence: DNA sequence (str)
    :param k: length of the kmers (int)

    :return: tuple containing the kmer frequency dictionary, the reassembled
    sequence, and the percentage of similarity between the original sequence
    and the reassembled sequence.
    """
    # call the kmer_breaks function to retrieve the kmer list
    kmer_dict = kmer_breaks(sequence, k)
    # create empty variables to be filled
    kmers_to_assemble = {}
    reassemble_list = []
    start = ''
    result_sequence = ''
    # for each kmer, find the overlapping kmers
    for key1 in kmer_dict:
        # kmers_to_assemble is composed by each kmer as a key, and a list with
        # the overlapping options from before and after the kmer (one base
        # before and one after).
        kmers_to_assemble[key1] = []
        for key2 in kmer_dict:
            if key1 != key2:
                if key2[1:] == key1[0: k - 1]:
                    kmers_to_assemble[key1].append(key2)
                if key1[1:] == key2[0: k - 1]:
                    kmers_to_assemble[key1].append(key2)
    # using the kmers_to_assemble dictionary, we join the kmer with the
    # overlapping kmers.
    for key in kmers_to_assemble:
        # if there is only one kmer overlapping this kmer is the start or the
        # end of the sequence.
        if len(kmers_to_assemble[key]) < 2:
            value = kmers_to_assemble[key][0]
            # check if the overlapping is in the end of the kmer, if yes we
            # found the starting kmer.
            if key[1:] == value[0:k - 1]:
                start += key + value[-1]
        # if there are two overlapping possibilities, they are the beginning
        # and the end of the key kmer. Check which one is which and join them.
        elif len(kmers_to_assemble[key]) == 2:
            value1 = kmers_to_assemble[key][0]
            value2 = kmers_to_assemble[key][1]
            if key[1:] == value1[0:k - 1]:
                reassemble_list.append(value2[0] + key[0:] + value1[-1])
            else:
                reassemble_list.append(value1[0] + key[0:] + value2[-1])
        # if there are more then 2 options, there are repetition, so save all
        # the overlapping possibilities.
        else:  # if len(kmers_to_assemble) > 2
            for kmer in kmers_to_assemble[key]:
                value = kmer
                if key[1:] == value[0:k - 1]:
                    reassemble_list.append(key + value[-1])
                else:
                    reassemble_list.append(value[0] + key)
    # if the starting kmer was found, use it as the first kmer in the result
    # sequence
    if start:
        result_sequence += start
        index = 0
    # if not, choose a random kmer to start.
    else:
        result_sequence += random.choice(reassemble_list)
        index = 1

    while len(result_sequence) < len(sequence):
        # while the result sequence is smaller than the original, keeps looping
        # through the kmers.
        for kmer in reassemble_list:
            # when an overlap is found, add the last base to the result
            # sequence.
            if result_sequence[index:] == kmer[0:len(kmer) - 1]:
                result_sequence += kmer[-1]
                index += 1
                break
        # if the loop ends without a match, add a random kmer (this avoids an
        # infinite loop)
        else:
            result_sequence += reassemble_list[index]
            # add k to index, so the overlapping check is done in the last
            # kmer of the sequence.
            index += k
        # if index is higher than the size of the sequence, break the while loop
        # (also to avoid infinite loop)
        if index > len(sequence):
            break
    # make sure that both sequences have the same length.
    result_sequence = result_sequence[0:len(sequence)]

    # calculate the similarity, by comparing each bases.
    matchs = 0
    sequence = sequence.upper()
    for i in range(len(sequence)):
        if sequence[i] == result_sequence[i]:
            matchs += 1
    percentage = round(matchs / len(sequence) * 100, 2)

    return kmer_dict, result_sequence, percentage


def visualize_graph(nodes, output):
    """
    Function to plot the De Bruijn graph in a pdf file.

    :param nodes: dictionary containing the possible kmers
    :param output: output pdf file name.

    :return: create the pdf file containing the The Bruijn Graph.
    """
    # remove '.pdf' from the name, it is included automatically
    output = output.strip('.png')
    # create the graph variable
    db_graph = Digraph(comment='De Bruijn Graph', format='png')
    # find all the kmers that are overlapping and save it in a dictionary
    kmers_to_assemble = {}
    for key1 in nodes:
        kmers_to_assemble[key1] = []
        for key2 in nodes:
            if key1 != key2:
                if key2[1:] == key1[0: len(key1) - 1]:
                    kmers_to_assemble[key1].append(key2)
                if key1[1:] == key2[0: len(key2) - 1]:
                    kmers_to_assemble[key1].append(key2)
    # for each node add the arrows
    for node in kmers_to_assemble:
        for kmer in kmers_to_assemble[node]:
            value = kmer
            if node[1:] == value[0:len(value) - 1]:
                # nodes are (k-1)-mers
                combination = node[0:]
                db_graph.edge(node[:len(node)-1], value[:len(value)-1],
                              label=combination, fontcolor='red')

    db_graph.render(output, view=True, cleanup=True)


def retrieve_sequence(file_path):
    """
    Function to retrieve a sequence in fasta format from a file. Only the first
    sequence will be retrieved if the file contains more than one.
    :param:
        file_path (str): Path of the fasta file.

    :return:
        String containing only the sequence.
    """
    # First checks if the file exists
    if os.path.exists(file_path):
        # Open the file, search for a line starting with ">" and return the next
        # line.
        sequence = ''
        with open(file_path, 'r') as file:
            if file.readline().startswith('>'):
                for line in file:
                    if not line.startswith('>'):
                        sequence += line.strip()
                    else:
                        break
        return sequence
    else:  # if the file does not exist.
        print('The file was not found.')


################################################################################
# Tkinter gui functions
# These functions are only for the GUI functionality.


def retrieve_input():
    """
    Function to retrieve the sequence and kmer input from the text boxes.
    The function result_window is called.
    """
    # Check if both information are passed (sequence and kmer)
    if sequence_box.get("1.0", 'end-1c') and kmer_selection.get("1.0",
                                                                'end-1c'):
        # retrieve the sequence and kmer.
        sequence = sequence_box.get("1.0", 'end-1c').upper()
        kmer = kmer_selection.get("1.0", 'end-1c').strip()
        # check if kmer is a number.
        if kmer.isnumeric():
            # call the result_window function.
            result_window(sequence, int(kmer))


def open_file():
    """
    Function for the Open file button, that opens a file dialog to let the user
    select a file.
    The function calls the retrieve_sequence function, and the returned sequence
    is inserted in the text box.
    """
    # Open the file dialog window.
    file = filedialog.askopenfilename()
    # If file is selected. If user clicks cancel, nothing happens.
    if file:
        # Call the retrieve_sequence function with the selected file.
        sequence = retrieve_sequence(file)
        if sequence:  # if the function returns a sequence, insert in the text
            # box.
            sequence_box.insert(tk.END, sequence)
        else:  # if the function returns False, a warning message appears in
            # the text box
            sequence_box.insert(tk.END, 'File not in fasta format.')


def result_window(sequence, kmer):
    """
    Function for the Submit button, that destroys the initial window and creates
    the result window, where the result information will be displayed.
    In this function the following are created:
        Main label: holds everything to facilitate destruction.
        Kmer box: List box displaying all the possible kmer and their frequency.
        Kmer box scrollbars: both x and y axis scrollbars.
        Result box: Text box displaying the original sequence and the result
            sequence for comparison.
        Result box scrollbar: y axis scrollbar for the result box.
        Similarity box: Box displaying the similarity percetage between the
            original sequence and the result sequence.
        New sequence button: allow a new sequence input.
        Close button: close the program.

    :param sequence: DNA sequence (str)
    :param kmer: length of the kmers (int)

    The function calls the reassembler function.
    If the draw_graph checkbox is selected, it opens a save file window, and
    calls the visualize_graph function.
    """
    global main_label
    nodes, result, similarity = reassembler(sequence, kmer)
    similarity = str(similarity) + '%'

    # Destroy the actual main label
    main_label.destroy()

    # Recreate the main label
    main_label = tk.Label(main_window, background='gray30')
    main_label.grid(row=0, column=0, padx=30, pady=10)

    # Kmer list box and scroll bar
    kmer_label = tk.Label(main_label, background='gray30')
    kmer_label.grid(row=0, column=0, rowspan=5, padx=10)
    kmer_box_text = tk.Label(kmer_label, text='Kmer List', background='gray60')
    kmer_box_text.grid(row=0, column=0, sticky='w')
    kmer_box = tk.Listbox(kmer_label, height=25, width=45)
    kmer_box.grid(row=1, column=0, columnspan=3)

    # y axis scrollbar for the kmer box
    kmer_box_scrollv = tk.Scrollbar(kmer_label, orient=tk.VERTICAL,
                                    command=kmer_box.yview)
    kmer_box_scrollv.grid(row=1, column=4, sticky='nsw')
    kmer_box.config(yscrollcommand=kmer_box_scrollv.set)

    # x axis scrollbar for the kmer box
    kmer_box_scrollh = tk.Scrollbar(kmer_label, orient=tk.HORIZONTAL,
                                    command=kmer_box.xview)
    kmer_box_scrollh.grid(row=6, column=0, columnspan=3, sticky='new')
    kmer_box.config(xscrollcommand=kmer_box_scrollh.set)

    # Populate the kmer box
    for node in nodes:
        if nodes[node] > 1:
            kmer_box.insert(tk.END, node + ' (x' + str(nodes[node]) + ')\n')
        else:
            kmer_box.insert(tk.END, node + '\n')

    # Create text box for the result and display the result.
    result_label = tk.Label(main_label, text='Reassembled sequence',
                            background='gray60')
    result_label.grid(row=0, column=1, sticky='ws')
    result_box = tk.Text(main_label, height=19, width=100)
    result_box.grid(row=1, column=1, columnspan=5, sticky='nwes')
    result_box.insert(tk.END, 'Original sequence:\n' + sequence)
    result_box.insert(tk.END, '\n\nResult sequence:\n' + result)

    # Result box scrollbar
    result_box_scroll = tk.Scrollbar(main_label, orient=tk.VERTICAL,
                                     command=result_box.yview)
    result_box_scroll.grid(row=1, column=6, sticky='nsw')
    result_box.config(yscrollcommand=result_box_scroll.set)

    # Similarity box
    similarity_label = tk.Label(main_label, background='gray60')
    similarity_label.grid(row=2, column=1, sticky='nw', pady=5)
    similarity_text = tk.Label(similarity_label,
                               text='Similarity between sequences',
                               background='gray60')
    similarity_text.grid(row=0, column=0, sticky='nws')
    similarity_result = tk.Label(similarity_label, text=similarity,
                                 background='white')
    similarity_result.grid(row=1, column=0, sticky='nwse')

    # New sequence button
    new_sequence_button = tk.Button(main_label, text='New Sequence',
                                    command=restart_window,
                                    background='gray60')
    new_sequence_button.grid(row=3, column=5, sticky='ew')

    # Close button
    close_button = tk.Button(main_label, text='Close',
                             command=main_window.destroy,
                             background='gray60')
    close_button.grid(row=4, column=5, sticky='we')

    # if checkbox is selected
    if var1.get():
        # open a save file window and call the visualize_graph function.
        output = filedialog.asksaveasfilename()
        visualize_graph(nodes, output)


def clear_box():
    """
    Function for the clear button, that clears the text box.
    """
    sequence_box.delete('1.0', tk.END)
    kmer_selection.delete('1.0', tk.END)


def initial_window():
    """
    Function that creates everything except the text boxes in the initial
    window.
    In this function the following are created:
        Scrollbar for the sequence box: y axis scrollbar
        Open file button: calls the open_file function.
        Submit button: call result_window function.
        Draw graph check box: make draw graph True.
        Clear button: call clear_box function.
        Close button: close the program.
    """
    # Text information over the text box.
    insert_label = tk.Label(main_label, text='Insert sequence in the box',
                            background='gray60')
    insert_label.grid(row=0, column=0, sticky='w')

    # Scrollbar for sequence box
    seq_box_scroll = tk.Scrollbar(main_label, orient=tk.VERTICAL,
                                  command=sequence_box.yview)
    seq_box_scroll.grid(row=1, column=5, sticky='nsw')
    sequence_box.config(yscrollcommand=seq_box_scroll.set)

    # Open file Button
    open_file_button = tk.Button(main_label, width=7, text='Open file',
                                 command=open_file, background='gray60')
    open_file_button.grid(row=2, column=3, sticky='we')

    # Submit button
    submit_button = tk.Button(main_label, width=7, text='Submit',
                              command=retrieve_input, background='gray60')
    submit_button.grid(row=2, column=4, sticky='we')

    # Clear button
    clear_button = tk.Button(main_label, text='Clear',
                             command=clear_box, background='gray60')
    clear_button.grid(row=4, column=3, sticky='we')

    # Close button
    close_button = tk.Button(main_label, text='Close',
                             command=main_window.destroy, background='gray60')
    close_button.grid(row=4, column=4, sticky='we')

    # Draw graph checkbox
    check_box = tk.Checkbutton(main_label, text='Draw Graph',
                               background='gray60', variable=var1, onvalue=1,
                               offvalue=0)
    check_box.grid(row=2, column=2, sticky='w', padx=10)


def restart_window():
    """
    Function for the New sequence button, that destroys the main label and
    create the initial window again.
    In this function, the following are created:
        Main Label.
        Sequence box.
        Kmer box
    Calls the initial_window function to create the other features.
    """
    # Use global variables to avoid error when retrieving the information in
    # the input boxes.
    global main_label
    global sequence_box
    global kmer_selection

    # destroy the main label.
    main_label.destroy()

    # recreate the main label
    main_label = tk.Label(main_window, background='gray30')
    main_label.grid(row=0, column=0, padx=30, pady=10)

    # Textbox
    sequence_box = tk.Text(main_label, height=15, width=140)
    sequence_box.grid(row=1, column=0, columnspan=5, sticky='nwes')

    # Kmer selection box
    kmer_label = tk.Label(main_label, background='gray60')
    kmer_label.grid(row=2, column=0, sticky='nw', pady=5)
    kmer_frame = tk.Label(kmer_label, text='Kmer Size',
                          background='gray60')
    kmer_frame.grid(row=0, column=0, sticky='nw')
    kmer_selection = tk.Text(kmer_label, height=1, width=7)
    kmer_selection.grid(row=1, column=0, sticky='nw')

    # call the initial_window function for the other features.
    initial_window()


################################################################################
# Running code

# If -c is used, run the code with the information in argparse.
if args.command:
    # check if all mandatory information are passed.
    if args.input_sequence or args.infile and args.kmer:
        k = args.kmer
        graph = args.draw_graph
        # check if sequence is passed manually or a file.
        if args.input_sequence:  # sequence
            dna_sequence = args.input_sequence.upper()
            kmers, result_seq, similarity_per = reassembler(dna_sequence, k)
        elif args.infile and args.kmer:  # file
            dna_sequence = retrieve_sequence(args.infile)
            if dna_sequence:  # if retrieve_sequence return a sequence
                kmers, result_seq, similarity_per = reassembler(dna_sequence, k)
            else:  # if retrieve_sequence returns False
                print('File not in fasta format.')

        # Print the results.

        # Kmer list and the frequencies
        print('\nKmer List:')
        for kmer in kmers:
            if kmers[kmer] > 1:
                print(kmer + ' (x' + str(kmers[kmer]) + ')')
            else:
                print(kmer)
        # Original sequence and result sequence for comparison
        print('\nOriginal sequence:\n' + dna_sequence)
        print('\nReassembled sequence:\n' + result_seq)
        # Similarity percentage between the sequences
        print('\nThe similarity between the original sequence and the result \
sequence is ' + str(similarity_per) + '%')

        # if --draw_graph flag is used, call the function to create the pdf.
        if args.draw_graph:
            if args.outfile:  # if output name is given
                output_file = args.outfile
            else:  # output default name
                output_file = 'result_de_bruijn_graph.png'
            visualize_graph(kmers, output_file)

    else:  # Help message if not all needed information is found.
        print('\nPlease submit a sequence or a file, and the kmer size.\n\
Use -h for more information.')

# If -c not used, create the interface
if not args.command:
    # Create the main window of the GUI.
    main_window = tk.Tk()
    main_window.title('De Bruijn Assembler')
    main_window.geometry('1200x480')
    main_window.configure(background='gray30')

    # Create the first text boxes.
    main_label = tk.Label(main_window, background='gray30')
    main_label.grid(row=0, column=0, padx=30, pady=10)

    # Textbox
    sequence_box = tk.Text(main_label, height=15, width=140)
    sequence_box.grid(row=1, column=0, columnspan=5, sticky='nwes')

    # Kmer selection box
    kmer_label = tk.Label(main_label, background='gray60')
    kmer_label.grid(row=2, column=0, sticky='nw', pady=5)
    kmer_frame = tk.Label(kmer_label, text='Kmer Size',
                          background='gray60')
    kmer_frame.grid(row=0, column=0, sticky='nw')
    kmer_selection = tk.Text(kmer_label, height=1, width=7)
    kmer_selection.grid(row=1, column=0, sticky='nw')

    var1 = tk.IntVar()

    # Call the function to complete the initial window.
    initial_window()

    # Start the main loop.
    main_window.mainloop()
