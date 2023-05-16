#!/bin/env/python3
import itertools
import numpy as np
import re
import math
import sys
from molmass import Formula
ionized_mass = {
    '[M+H]+': 1.00727646677,
    '[M-H]-': -1.00727646677,
    '[M+Na]+':22.989218,
    '[M+NH4]+':18.033823
}

atom_mass = {
    "H": 1.007825,
    "C": 12.00000,
    "N": 14.003074,
    "O": 15.994915,
    'F': 18.99840322,
    'Na':22.989769282,
    "P": 30.973762,
    "S": 31.972071,
    "Cl": 34.968852721,  # Cl 37: 36.9659
    "Br": 78.9183,  # Br 81: 80.9163
    'I': 126.9045

    # 'Co': 58.93319429
}
atom_list = ['H', 'C', 'N', 'O', 'F', 'Na','P', 'S', 'Cl', 'Br', 'I']
numpy_formula_format = np.int16

atom_dict = {a: i for i, a in enumerate(atom_list)}
len_atom_dict = len(atom_dict)
atom_mass_array = np.zeros(len_atom_dict, np.float32)
for atom in atom_dict:
    atom_mass_array[atom_dict[atom]] = atom_mass[atom]


class MolecularFormula(object):
    __slots__ = ['_data', '_hash']

    def __init__(self, data=None):
        if data is not None:
            self._data = np.array(
                data, numpy_formula_format, copy=True)
        else:
            self._data = np.zeros(len_atom_dict, numpy_formula_format)
        self._hash = None
        pass

    def __getitem__(self, item):
        return self._data[atom_dict[item]]

    def __setitem__(self, key, value):
        self._data[atom_dict[key]] = value

    def __str__(self):
        string = ''

        for atom in atom_list:
            atom_num = self[atom]
            if atom_num:
                if atom_num > 1:
                    string += atom + str(atom_num)
                else:
                    string += atom
        return string

    def from_string(self, key_string):
        f = Formula(key_string)
        # f.composition()
        all_atom_nums = []
        for i in range(len(f.composition())):
            all_atom_nums.append((f.composition()[i][0], f.composition()[i][1]))
        for atom_num in all_atom_nums:
            self[atom_num[0]] = int(atom_num[1])
        pass

    def get_data(self):
        return self._data

    def get_mass(self):
        return np.sum(atom_mass_array * self._data)


def mass_to_formula(mass, mass_error, addition, precursor_formula=None):
    mass = float(mass)
    mass_error = float(mass_error)
    if precursor_formula is None:
        return precursor_mass_to_formula(mass, mass_error, addition)
    else:
        mol = MolecularFormula()
        mol.from_string(precursor_formula)
        result = product_mass_to_formula(mass, mass_error, addition, mol)
        return result

def nl_to_formula(mass, mass_error, molecular_formula):
    mass = float(mass)

    mass_error = float(mass_error)

    mol = MolecularFormula()
    mol.from_string(Formula(molecular_formula).formula)
    # if Formula(molecular_formula).isotope.mass <=200:
    #     mass_error = 0.01
    # if Formula(molecular_formula).isotope.mass >200:
    #     mass_error = Formula(molecular_formula).isotope.mass*mass_error/1E6
    result = nl_mass_to_formula(mass, mass_error, mol)
    formula = []
    formula_mass = []
    for r in result:
        formula.append(r.__str__())
        # formula_mass.append(r.get_mass())
    return(formula)
        # return(r.__str__())
        # break
        # return(r)
        # break

def nl_mass_to_formula(mass, mass_error, precursor_formula):

    lo_mass = mass - mass_error
    hi_mass = mass + mass_error
    precursor_data = precursor_formula.get_data()
    formula_range = [range(x + 1) for x in precursor_data]
    all_possible_candidate_formula = np.array(
        list(itertools.product(*formula_range)), numpy_formula_format)
    all_possible_mass = np.sum(
        atom_mass_array * all_possible_candidate_formula, axis=1)
    candidate_data = all_possible_candidate_formula[(lo_mass <= all_possible_mass) & (all_possible_mass <= hi_mass)]
    result = []
    for data in candidate_data:
        formula = MolecularFormula(data)
        result.append(formula)
    return result
# def _calculate_formula(mass_start, mass_end, candidate_formula_array, cur_i, result):
#     atom_mass_cur = atom_mass_array[cur_i]
#     atom_num = math.floor(mass_end / atom_mass_cur)
#     if cur_i == 0:
#         # This is H
#         h_num_low = mass_start / atom_mass_cur
#         if atom_num >= h_num_low:
#             candidate_formula_array[0] = atom_num
#             result.append(MolecularFormula(candidate_formula_array))
#     else:
#         for i in range(atom_num):
#             f = np.copy(candidate_formula_array)
#             f[cur_i] = i
#             _calculate_formula(mass_start - i * atom_mass_cur, mass_end - i * atom_mass_cur,
#                                f, cur_i - 1, result)
# if __name__ == "__main__":
#     if len(sys.argv) == 4 or len(sys.argv) == 5:
#         result = mass_to_formula(*sys.argv[1:])
#         for r in result:
#             print(r, r.get_mass())
#     else:
#         print("""
# python3 mass_to_formula.py ion_mass mass_error ion_type [precursor_formula]
#
# For example:
#
# python3 mass_to_formula.py 508.003 0.002 [M+H]+
#
# or
#
# python3 mass_to_formula.py 136.0616 0.002 [M+H]+ C10H16N5O13P3
#
# """)
