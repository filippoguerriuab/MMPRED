# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:19:09 2022

@author: pippo
"""

import pandas as pd
import numpy as np
import sys


def compute_ranks_from_rand_pred(in_file, dict_filename):


    df_rnd_scores = pd.read_csv(in_file)
    hla_names = df_rnd_scores.columns.to_list()
    hla_names = hla_names[4::]
    ranks_dict = {}
    for hla in hla_names:
        scores = np.array([round(float(x),3) for x in df_rnd_scores[hla] if x != None])
        ranks_dict[hla] = scores


    out_cont = ''
    for hla in ranks_dict:
        out_cont += hla.replace('pred_outcome_HLA_', "") + "," + ",".join([str(x) for x in ranks_dict[hla]]) + "\n"

    f = open(dict_filename, 'w')
    f.write(out_cont)


def apply_ranks(ranks_dict, hla, score):

    score_type = type(score)
    print(score_type)
    if score_type == int or score_type == float:
        rank = np.sum(ranks_dict[hla]>score)/len(ranks_dict[hla])*100
        return rank

    elif score_type == list or score_type == np.ndarray:
        if score_type == list:
            score = np.array(score)
        rank = np.array(len(score))
        for i, s in enumerate(score):
            rank[i] = np.sum(ranks_dict[hla]>s)/len(ranks_dict[hla])*100



if __name__ == '__main__':
    in_file = sys.argv[1]
    dict_filename = sys.argv[2]
    compute_ranks_from_rand_pred(in_file, dict_filename)
