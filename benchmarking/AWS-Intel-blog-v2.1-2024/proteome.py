import sys
import os
import re
import numpy as np


with open("./uniprotkb_proteome.fasta", "r") as f:
    lines = f.readlines()
    i = -1
    protien_list = []
    proteome_list = []
    for line in lines:
        if ">sp|" in line:
            if i >= 0:
                proteome_list.append(protien_list)
                protien_list = []
            i = i + 1
            protien_list.append(line)
        else:
            protien_list.append(line)
    proteome_list.append(protien_list)
    
    sorted_list = sorted(proteome_list, key=lambda x: len(''.join(x[1:])), reverse=False)
    i = 0
    sum = 0
    small_db = open("short_db", "r")
    small_list = [line.rstrip() for line in small_db.readlines()]
    small_db.close()
    #lines = small_db.readlines()
    #print(lines)

    long_db = open("long_db", "r")
    long_list = [line.rstrip() for line in long_db.readlines()]
    long_db.close()
    os.mkdir("~/celegans_samples")
    os.mkdir("~/celegans_samples_long")
    for pl_list in sorted_list:
        sum = sum + len(''.join(pl_list[1:]))
        print(i, len(''.join(pl_list[1:])))
        if "celegans_"+str(i)+".fa" in small_list:
            with open("~/celegans_samples/celegans_" + str(i) + ".fa", "w") as f:
                f.writelines(pl_list)
        
        if "celegans_"+str(i)+".fa" in long_list:
            with open("~/celegans_samples_long/celegans_" + str(i) + ".fa", "w") as f:
                f.writelines(pl_list)



        i = i + 1
        
    print(sum/i)
