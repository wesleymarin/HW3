import numpy as np
import glob
import os
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import sys
import pandas as pd
import math
from multiprocessing import Pool
import pickle
from sklearn import metrics

def read_sequences(dir):
    """
    Creates a list of sequence objects from directory 'dir'

    Input: Directory where fasta sequences are held
    Output: List of sequence objects
    """
    ## Only get files that end in .fa
    files = glob.glob(dir + '/*.fa')

    ## Initialize sequnce object list
    sequences = []

    ## Send all found files to the read_sequence function
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):
        sequences.append(read_sequence(filepath))

    return(sequences)


def read_sequence(filepath):
    """
    Creates a sequence object from a file.

    Input: Filepath to create sequence from
    Output: Sequence object
    """

    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    ## Make sure it is a fasta file
    if name[1] != ".fa":
        raise IOError("%s is not a Fasta file"%filepath)

    ## Initialize the sequence object with the filename (without '.fa')
    sequence = Sequence(name[0])

    ## Start up residue numbering
    res_num = 1
    with open(filepath, "r") as f:
        for line in f:

            ## Strip newlines
            line = line.rstrip()

            ## Dont read in the identity line
            if line [0] != '>':

                ## For character in the line, add the char to sequence.residues, along with
                ## the residue number
                for i in range(0, len(line)):
                    sequence.residues.append(Residue(line[i], res_num))
                    res_num += 1

    return(sequence)


class Sequence:
    def __init__(self, name):
        self.name = name
        self.residues = []

    def __repr__(self):
        return self.name

class Residue:

    def __init__(self, type, number):
        self.type = type
        self.number = number

    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)


def read_match_matrix(dir, name):
    """
    Creates a pandas dataframe from a scoring matrix.

    Input: Directory where the scoring matrix is found
           Name of the scoring matrix file

    Output: Pandas dataframe with scores
    """

    ## Find the file
    file = glob.glob(dir + "/" + name)

    ## Pointless loop that only goes once
    for filepath in file:

        with open(filepath, 'r') as f:

            ## Need to be able to tell when the first line happens to get the index names
            line_count = 1

            for line in f:
                ## Take out the newlines and split every line into single elements
                line = line.rstrip()
                line = line.split()

                ## If it is not a comment line
                if line[0] != '#':

                    ## Use the first line to initialize the dataframe with index and column names
                    if line_count == 1:
                        match_matrix = pd.DataFrame(np.zeros((len(line), len(line)), np.int8), index=line, columns=line)
                    else:

                        ## Throw all of the values into the dataframe, cast values as ints
                        match_matrix.iloc[line_count - 2] = [int(l) for l in line]

                    line_count += 1

    return match_matrix

def read_match_pairs(dir, name):
    """
    Creates a list of pairs, which are themselvs lists. Pair names come from a file.

    Input: Directory where pair file is found
           Name of pair file
    Output: List of lists
    """
    ## Get the file path, this time I took out the pointless for loop
    file = glob.glob(dir + "/" + name + ".txt")[0]

    with open(file, 'r') as f:
        ## Initialize the list of pairs
        pair_list = []

        for line in f:

            ## Take out the new lines, then split each line into two using " " as a delim
            line = line.rstrip()
            line = line.split(" ")

            ## Initialize the pair list
            pair = []

            ## Each line is a pair
            for element in line:

                ## Take out the needless parts of the pair names
                element = element.split("/")[1]
                element = element.split(".")[0]

                ## Add each name to the pair
                pair.append(element)

            ## Add the pair to the list of pairs
            pair_list.append(pair)
    return(pair_list)

def match_score(match1, match2, match_matrix):
    """
    Find the score of matching match1 with match2 in the match_matrix

    Input: Index name to match (residue type)
           Column name to match (residue type)
           Dataframe with scores
    Output: Score (int)
    """
    try:
        return(match_matrix[match1.upper()][match2.upper()])
    except KeyError:
        print("Whoops! Those residues were not found in the match matrix.")
        print(match1)
        print(match2)


def create_score_table(sequence1, sequence2, match_matrix, gap_open, gap_cont):
    """
    Based on the residue match scores found in match_matrix, this function computes a score table for
    every single residue combination. These scores are computed from the top left of the table to the bottom right.

    Along with the scores, a direction vector is included in every cell, indicating the path of best score propogation.

    Input: First sequence object to match
           Second sequence object to match
           DataFrame of residue match scores
           Gap opening penalty
           Gap continue penalty

    Output: DataFrame with len(sequence1)+1 rows and len(sequence2)+1 columns, with scores for each combination.
    """

    ## Initialize a general gap_score variable that we can change for opening and continuing
    gap_score = gap_open

    ## Initialize a datatype for our dataframe cells for an int score and a direction vector
    dt = np.dtype([('score', np.int8), ('direction', str, 3)])

    ## Initialize the np array to build the dataframe from, extra row and column are for pure gap scores
    initial_table = np.empty((len(sequence1.residues)+1, len(sequence2.residues)+1), dtype=dt)

    ## Initialize the score dataframe with sequence.residue numbers as index and col names
    score_table = pd.DataFrame(initial_table,
                               index=["gap"] + [resi.number for resi in sequence1.residues],
                              columns=["gap"] + [resi.number for resi in sequence2.residues],
                              dtype=None)

    ## Initialize the top right score to be 0 with no direction vector
    score_table.iloc[0][0] = (0, '---')

    ## Initialize the gap row and column, which would represent an alignment of pure gaps
    for residue1 in sequence1.residues:

        ## For the gap column, we want to make all the direction vectors point up
        if residue1.number == 1:

            ## If it is the first gap score, use gap_open
            score_table.iloc[residue1.number][0] = (gap_open, '-u-')
        else:

            ## If not the first gap score, add the previous row score to the gap continue score
            score_table.iloc[residue1.number][0] = (score_table.iloc[residue1.number-1][0][0] + gap_cont, '-u-')


    ## Do the same things for the row of pure gaps
    for residue2 in sequence2.residues:
        if residue2.number == 1:
            score_table.iloc[0][residue2.number] = (gap_open, '--l')
        else:
            score_table.iloc[0][residue2.number] = (score_table.iloc[0][residue2.number-1][0] + gap_cont, '--l')



    ## This is the meat of the algorithm, iterate over every combination of residues in sequence1 and sequence2
    for residue1 in sequence1.residues:
        for residue2 in sequence2.residues:

            ## Initialize a score list with the score of the previous diagonal match + the score of the new match
            ## The direction vector for this option would be diagonal because that was the movement it represents
            score_list = [(score_table.iloc[residue1.number-1][residue2.number-1][0] + match_score(residue1.type,
                                                                                               residue2.type,
                                                                                               match_matrix), 'd--')]


            ## Before computing the up direction gap scores,
            ## check to see if there already is a gap in the up direction
            if 'u' in score_table.iloc[residue1.number-1, residue2.number][1]:
                gap_up = gap_cont
            else:
                gap_up = gap_open

            ## Add a score for moving down (with a vector that points back up)
            score_list.append((score_table.iloc[residue1.number-1, residue2.number][0] + gap_up, '-u-'))

            ## Do the same thing for the left direction gap scores
            if 'l' in score_table.iloc[residue1.number, residue2.number-1][1]:
                gap_left = gap_cont
            else:
                gap_left = gap_open

            ## Add a score for moving right (with a vector that points back left)
            score_list.append((score_table.iloc[residue1.number, residue2.number-1][0] + gap_left, '--l'))

            ## Add a zero score with no direction vector (for smith-waterman algorithm)
            score_list.append((0, '---'))

            ## Out of all the scores in the score list, find the max
            greatest_score = max([score[0] for score in score_list])

            ## Grab all scores + directions for every score that equals the max score
            greatest_scores_raw = [score for score in score_list if score[0] == greatest_score]

            ## Initialize a final direction string
            direction_string = list('---')

            ## For every greatest score found, add the direction associated with the score to the final
            ## direction string
            for score in greatest_scores_raw:
                if score[1] == 'd--':
                    direction_string[0] = 'd'
                elif score[1] == '-u-':
                    direction_string[1] = 'u'
                elif score[1] == '--l':
                    direction_string[2] = 'l'
                else:
                    None

            ## Turning the direction string back into a string
            direction_string = "".join(direction_string)

            ## Put the greatest score and direction string together in a tuple
            score_entry = (greatest_score, direction_string)

            ## Put the tuple into the score table for this residue combination
            score_table[residue2.number][residue1.number] = score_entry

    return(score_table)


def create_match_strings(score_table, sequence1, sequence2):
    """
    From a score table computed by create_score_table, create a match between sequence1 and sequence2 by finding
    the max score in the score table and following the direction vectors back to 0, saving all the matches along
    the way.

    Input: Dataframe of match scores and directions
           Sequence object to match
           Sequence object to match

    Output: list of matched residues for sequence1
            list of matched residues for sequence2
            Max score from which the matches were started from
    """

    ## Find the largest score in the table, take the first if there are multiple
    largest_score = score_table.max(axis=1).max()[0]

    ## Initialize a furthest match location (to find the match that is closest to the bottom right)
    furthest_match = (0, 0)

    ## Iterate through all of the rows and columns
    for row in score_table.index:
        for column in score_table.columns:

            ## If any score equals the largest score and has a further position than any other largest score,
            ## save the position of this score
            if score_table[column][row][0] == largest_score:
                if column + row > furthest_match[0] + furthest_match[1]:
                    furthest_match = (row, column)


    ## Initialize the residue match strings
    match_string1 = []
    match_string2 = []

    ## Initialize the starting match location to be the furthest match location
    match_coord = furthest_match

    ## Initialize the current score to be the largest score
    current_score = largest_score

    ## Follow the match directions until the score falls to 0
    while current_score > 0:

        ## Get the value of the score table cell for the match coordinates
        current_match = score_table[match_coord[1]][match_coord[0]]

        ## Extract the direction vector as a list of characters
        current_direction = list(current_match[1])

        ## Replace current score with this new match score
        current_score = current_match[0]

        ## If there is a diagonal in the direction vector
        if current_direction[0] == 'd':

            ## Add the residues at the current coordinates to the match strings
            match_string1.append(sequence1.residues[match_coord[0]-1].type)
            match_string2.append(sequence2.residues[match_coord[1]-1].type)

            ## Update the coordinate to correspond with diagonal movement
            match_coord = (match_coord[0]-1, match_coord[1]-1)

        ## If there is an up in the direction vector
        elif current_direction[1] == 'u':

            ## Only add a matched residue to sequence1, the other match string gets a gap
            match_string1.append(sequence1.residues[match_coord[0]-1].type)
            match_string2.append("-")

            ## Move the coordinates up one
            match_coord = (match_coord[0]-1, match_coord[1])

        ## If there is a left in the direction vector
        elif current_direction[2] == 'l':

            ## Do the opposite of the last one
            match_string1.append("-")
            match_string2.append(sequence2.residues[match_coord[1]-1].type)

            match_coord = (match_coord[0], match_coord[1]-1)

        ## If there is no direction in the direction vector, pop the last element off because we went too far
        else:
            match_string1.pop()
            match_string2.pop()

        ## Bug fix where match coord got down to 0 and the index didnt work
        if match_coord[0] == 0:
            match_coord = ("gap", match_coord[1])
        if match_coord[1] == 0:
            match_coord = (match_coord[0], "gap")

    return (match_string1, match_string2, largest_score)
