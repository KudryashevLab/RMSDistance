#!/usr/bin/env python3
"""
fetch_changechains_align_save.py
Batch-fetch PDB models, rename chains to reference nomenclature,
align to reference structure, and save.
Usage: chimerax --nogui --script fetch_changechains_align_save.py > fetch_changechains_align_save.log 2>&1
"""
# -*- coding: utf-8 -*-

# ─── CONFIG ───────────────────────────────────────────────────────────────────

INPUT_FILE       = "models_with_state_short.txt"   # .txt list of PDB IDs to process
OUTPUT_DIR       = "./models"                       # where .cif files are saved
REFERENCE_MODEL  = "./reference_model.cif"         # reference model for alignment
REFERENCE_ID     = 1337                            # ChimeraX model ID for reference

# Alignment residue range (Cα atoms), TMD in case of RyR1
ALIGN_START = 4820
ALIGN_END   = 4920

# Sequence fingerprints for chain identification
SEQ_RYR = "QFLRTDDEVV"   # RyR1          (mandatory)
SEQ_CLS = "TISPGDGRTF"   # FKBP          (optional)
SEQ_CAM = "TEEQIAEFKE"   # Calmodulin    (optional)
SEQ_NAB = "QVQLQESGGG"   # Nanobody      (optional)

# ─── END CONFIG ───────────────────────────────────────────────────────────────

import os
import string
from chimerax.atomic import all_structures
from chimerax.core.commands import run

os.makedirs(OUTPUT_DIR, exist_ok=True)

abc = list(string.ascii_uppercase)

files = {}; cnt = 0
f = open(INPUT_FILE, 'r')
for i in f:
    files[cnt] = i.split()
    cnt = cnt + 1
f.close()

for i in files:
    if os.path.exists(OUTPUT_DIR + "/" + ''.join(files[i][0]) + ".cif"):
        continue
    run(session, "open " + ''.join(files[i][0] + " id " + str(i+1)))
    structures = all_structures(session)
    residues = structures[0].residues
    run(session, "sel sequence " + SEQ_RYR)
    ids = list(residues[residues.selected].unique_chain_ids)

    run(session, "sel sequence " + SEQ_CLS); cls = []
    if len(residues[residues.selected]) > 0:
        cls = list(residues[residues.selected].unique_chain_ids)
        ids = ids + cls

    run(session, "sel sequence " + SEQ_CAM); cam = []
    if len(residues[residues.selected]) > 0:
        cam = list(residues[residues.selected].unique_chain_ids)
        ids = ids + cam

    run(session, "sel sequence " + SEQ_NAB); nab = []
    if len(residues[residues.selected]) > 0:
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

    if len(cls) > 0:
        idx = idx + 4
        cls_new = []
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            cls_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids)) + idx - 4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(cls_new) + " " + ",".join(abc[idx-4:idx]))

    if len(cam) > 0:
        idx = idx + 4
        cam_new = []
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            cam_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids)) + idx - 4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(cam_new) + " " + ",".join(abc[idx-4:idx]))

    if len(nab) > 0:
        idx = idx + 4
        nab_new = []
        for j in abc[idx-4:idx]:
            run(session, "sel /" + j)
            run(session, "sel zone sel 10 residues 1 extend 0")
            nab_new.append(abc[abc.index(''.join(residues[residues.selected].unique_chain_ids)) + idx - 4])
        run(session, "changechains #" + str(i+1) + " " + ",".join(nab_new) + " " + ",".join(abc[idx-4:idx]))

    run(session, "open " + REFERENCE_MODEL + " id " + str(REFERENCE_ID))
    run(session, "align #" + str(i+1) + ":" + str(ALIGN_START) + "-" + str(ALIGN_END) + "@ca to #" + str(REFERENCE_ID) + ":" + str(ALIGN_START) + "-" + str(ALIGN_END) + "@ca")
    run(session, "save " + OUTPUT_DIR + "/" + ''.join(files[i][0]) + ".cif #" + str(i+1))
    run(session, "close #" + str(i+1))
    run(session, "close #" + str(REFERENCE_ID))
run(session, "quit")
