import pandas as pd
from molmass import Formula
from functools import reduce
import os
import numpy as np
from fuzzywuzzy import fuzz
def save_value_counts(data, column):
    vc = data[column].value_counts().rename_axis('unique_values').to_frame('counts')
    vc.index.name = 'unique_values'
    vc.reset_index(inplace=True)
    return(vc)
def check_missing_compound(col1, col2):
    missing = []
    if len(set(col1))<=len(set(col2)):
        shorter = set(col1)
        longer = set(col2)
    else:
        shorter = set(col2)
        longer = set(col1)
    for item in longer:
        if item in shorter:
            pass
        else:
            missing.append(item)
    return(missing)
def find_files(base, pattern):
    # print('i am in new new')
    score = []
    for filename in os.listdir(base):
        score.append(fuzz.ratio(filename,pattern))
    return(os.listdir(base)[np.argmax(score)])
def specify_column(keyword, column_name):
    score = []
    for name in column_name:
        score.append(fuzz.token_sort_ratio(keyword,name))
    return(column_name[np.argmax(score)])
def find_different_items(col1, col2):
    if len(col2)<len(col1):
        list_1 = col2
        list_2 = col1
    else:
        list_1= col1
        list_2 = col2
    main_list = list(set(list_2) - set(list_1))

    return(main_list)

# def get_mono_mass(data, formula_column = "Formula"):
#     mass = []
#     for index, row in data.iterrows():
#         mass.append(Formula(row[formula_column]).isotope.mass)
#     data.insert(data.columns.get_loc(formula_column)+1, "Monoisotopic mass", mass)
#     return(data)
def get_unique_list(list1):
    list_set = set(list1)
    unique_list = (list(list_set))
    return(unique_list)
def find_common_items(columns):
    return (list(reduce(set.intersection, map(set, columns))))




