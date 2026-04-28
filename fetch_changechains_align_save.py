files = {};cnt = 0
f = open('models_with_state.txt', 'r')
for i in f:
    files[cnt] = i.split()
    cnt = cnt+1
f.close()

from chimerax.atomic import all_structures
from chimerax.core.commands import run
import string
import os
abc = list(string.ascii_uppercase)
for i in files:
    if os.path.exists("./models/" + ''.join(files[i][0]) + ".cif"):
        continue
    run(session, "open " + ''.join(files[i][0] + " id " + str(i+1)))
    structures = all_structures(session)
    residues = structures[0].residues
    run(session, "sel sequence QFLRTDDEVV")
    ids = list(residues[residues.selected].unique_chain_ids)

    run(session, "sel sequence TISPGDGRTF");cls=[]
    if len(residues[residues.selected])>0:
        cls = list(residues[residues.selected].unique_chain_ids)
        ids = ids + cls

    run(session, "sel sequence TEEQIAEFKE");cam=[]
    if len(residues[residues.selected])>0:
        cam = list(residues[residues.selected].unique_chain_ids)
        ids = ids + cam

    run(session, "sel sequence QVQLQESGGG");nab=[]
    if len(residues[residues.selected])>0:
        nab = list(residues[residues.selected].unique_chain_ids)
        ids = ids + nab

    run(session, "changechains #" + str(i+1) + " " + ",".join(ids) + " " + ",".join(abc[0:len(ids)]))

    id_new = ['A']
    idx = 4
    for l in range(3):
        run(session, "sel /" + id_new[l] + ":1-100")
        run(session, "sel zone sel 10 residues 1 extend 0")
        id_new.append(''.join(list(set(''.join(residues[residues.selected].unique_chain_ids)) ^ set(id_new[l]))))
    run(session, "changechains #" + str(i+1) + " " + ",".join(id_new) + " " + ",".join(abc[0:idx]))


    if len(cls)>0:
        idx = idx+4
        cls_new = []
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            cls_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids))+idx-4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(cls_new) + " " + ",".join(abc[idx-4:idx]))

    if len(cam)>0:
        idx=idx+4
        cam_new=[]
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            cam_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids))+idx-4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(cam_new) + " " + ",".join(abc[idx-4:idx]))

    if len(nab)>0:
        idx=idx+4
        nab_new=[]
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            nab_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids))+idx-4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(nab_new) + " " + ",".join(abc[idx-4:idx])) 

    run(session, "open ../j126_final_lig_is.cif id 1337")
    run(session, "align #" + str(i+1) + ":4820-4920@ca to #1337:4820-4920@ca")
    run(session, "save ./models/" + ''.join(files[i][0]) + ".cif #" + str(i+1))
    run(session, "close #" + str(i+1))
    run(session, "close #1337")
