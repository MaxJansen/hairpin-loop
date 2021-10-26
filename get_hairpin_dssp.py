#!/usr/bin/python3

"""
What columns do I want?
pdb_filename, the secondary string output, Does it match hairpin pattern?,

"""
# Import packages
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import argparse
import os
import pandas
import re

original_wd = os.getcwd()

def line_read_dssp(ding):
    """
    Takes a pdb filepath and returns name and DSSP secondary structure output
    """
    # Break up lines from file to get pdb and path
    pdb_filename = ding.split('/')[-1]
    pdb_filename = pdb_filename.replace("\n", "")
    pdb_id = pdb_filename[0:4]
    pdb_dir = ding.split('/')[:-1]
    pdb_dir = '/'.join([str(n) for n in pdb_dir])

    # Go to directory and load pdb to get DSSP
    os.chdir(pdb_dir)
    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename, dssp="/work/upcorreia/bin/sequence/dssp")
    full_dssp = list(dssp)
    sec_str_result = [i[2] for i in full_dssp]
    sec_str_result = ''.join([str(n) for n in sec_str_result])

    return [pdb_filename, sec_str_result]

def quick_checks(sec_str_result):
    # Does the hairpin have a dash at the start and end?
    if sec_str_result[0] == '-' and sec_str_result[-1] == '-':
        start_finish = ''
    else:
        start_finish = "Start and end not like hairpin"

    # Is it longer than 10 residues?
    if len(sec_str_result) > 10:
        length10 = ''
    else:
        length10 = "10 or shorter"
    return [start_finish, length10]

def summ_dssp_anno(sec_str_result):
    """
    Returns summary. Categorises residue secondary struture
    assignment from DSSP into three main groups: Helix (H), Loop (L)
    and Sheet (B).
    """
    summ_dict = {'H': 'H', 'I': 'H', 'B': 'B', 'E': 'B', 'G': 'L', 'T': 'L',
                  'S': 'L', '-': 'L'}
    summ_list = list(map(lambda x: summ_dict[x], list(sec_str_result)))
    summ_str = ''.join(map(str,summ_list))

    return summ_str

def assess_summ(summ_str):
    """
    Take DSSP type summary, does it follow expected pattern?
    What percentage is loop?
    """
    # Does is follow a helix, loop, helix pattern?
    if re.search("^(L+H+L+H+L+)$", summ_str):
        summ_type = "Exact hairpin pattern"
    elif re.search("(L+H+L+H+L+)", summ_str):
        summ_type = "Contains hairpin pattern"
    else:
        summ_type = "Not a hairpin pattern"

    # Calculate percentage loop, subtract 2 for start and end
    full_len = len(summ_str)
    loop_counter = summ_str.count('L')-2
    pct_loop = loop_counter/full_len * 100

    return [pct_loop, summ_type]

def parser():
    """Retrieves the arguments from the command line."""

    input_file = argparse.ArgumentParser(description='Computes secondary \
    structure of hairpins and tells you what percentage is loop, performs checks.')
    input_file.add_argument('-in', required=True, metavar='input_file',
                      dest='inputFile', help='[-ab] ')
    arguments = input_file.parse_args()  # takes the arguments
    return [arguments]

#==============================
# Main commands
#===============================

# Get command line arguments, allows running the script for
# a file containing list of pdb's
COMMAND_LINE = parser()
source = COMMAND_LINE[0].inputFile
F = open(source, "r")
X = F.readlines()

# Get df ready to export results
column_names = ["pdb_filename", "sec_str_result", "pct_loop", "summ_type",
                "start_finish", "length10"]
df = pd.DataFrame(columns = column_names)

# Keep counting while running
all_hairpins = len(X)
i = 0

for line in X:
    # Run functions for each pdb
    dssp_result = line_read_dssp(line)
    summ = summ_dssp_anno(dssp_result)
    assessment = assess_summ(summ)
    checks = quick_checks(dssp_result[1])

    # concatenate results from function and add to final df
    x = np.array([dssp_result, assessment, checks])
    x = list(x.flatten())
    df.loc[len(df.index)] = x

    # Counter while running
    print(str(i) + " out of " + str(all_hairpins) + " done")
    i+=1

os.chdir(original_wd)
df.to_csv(r'dssp_loop_result.txt')
