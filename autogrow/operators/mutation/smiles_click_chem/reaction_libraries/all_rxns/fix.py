import os
import glob
import sys

path ="/home/jacob/Documents/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/temp.json"
new_path="/home/jacob/Documents/autogrow4/autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/temp_new.json"

printout = ""

new_rxn_num = 37
original_rxn_num = 1
with open(path,'r') as f:
    for line in f.readlines():
        if '": {' in line or "reaction_name" in line:
            line = line.replace("{}_".format(original_rxn_num),"{}_".format(new_rxn_num))
        elif '"RXN_NUM": ' in line:
            line = line.split(": ")[0] + ": {}\n".format(new_rxn_num)
            new_rxn_num = new_rxn_num + 1
            original_rxn_num = original_rxn_num + 1
        
        printout = printout + line

with open(new_path,'w') as f:
    f.write(printout)