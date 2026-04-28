from chimerax.atomic import all_structures
from chimerax.core.commands import run
import numpy as np
import glob
import os

files = {};cnt = 0
f = open('models_with_state.txt', 'r')
for i in f:
    files[cnt] = i.split()
    cnt = cnt+1
f.close()

with open("./ef12_s2s3.star", 'w') as f:
    f.write("\ndata_models\n\nloop_\n_ModelId\t#1\n_ModelState\t#2\n_Distance1\t#3\n_Distance2\t#4\n_Calstabin\t#5\n")
    for i in files:
        run(session, "open models/" + files[i][0] + ".cif" + f" id {i+1:d}")
        distance1 = run(session, f"distance #{i+1:d}/A:4101@ca #{i+1:d}/D:4730@ca");
        distance2 = run(session, f"distance #{i+1:d}/A:4075@ca #{i+1:d}/D:4736@ca");
        sel = run(session, "sel sequence TISPGDGRTF");cls=len(sel.models);
        f.write(files[i][0]+ "\t" + files[i][1] + f"\t{distance1:6.2f}\t{distance2:6.2f}\t{cls:d}\n")
        run(session, "close all")

