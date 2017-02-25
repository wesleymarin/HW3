from alignment_testing
import os
import matplotlib.pyplot as plt
import math
import numpy as np


def test_stuff():
    sequence_list = read_sequences("sequences/")
    match_matrix = read_match_matrix(".", "BLOSUM50")
    neg_pairs = read_match_pairs('.', 'Negpairs')
    pos_pairs = read_match_pairs('.', 'Pospairs')

    assert type(sequence_list[0]) is Sequence
    assert match_matrix['A']['A'] == 5
    assert len(neg_pairs) > 0 and len(neg_pairs) == len(pos_pairs)

    score_table = create_score_table(sequence_list[0], sequence_list[1], match_matrix, -5, -2)
    match1, match2, score = create_match_strings(score_table, sequence_list[0], sequence_list[1])

    assert type(score) is np.int64 and score > 0
