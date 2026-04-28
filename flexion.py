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

with open("./flexions.star", 'w') as f:
    f.write("\ndata_models\n\nloop_\n_ModelId\t#1\n_ModelState\t#2\n_FlexionTilt\t#3\n_FlexionRot\t#4\n_Shell2Tilt\t#5\n_Calstabin\t#6\n_RY1&2_angle\t#7\n_RY1&2_distance\t#8\n_Distance1\t#9\n_Distance2\t#10\n_PoreDistance\t#11\n_RY12_BSol_distance\t#12\n")
    for i in files:
        run(session, "open models/" + files[i][0] + ".cif" + f" id {i+1:d}")
        center = run(session, f"define centroid #{i+1:d}");
        axis_Z = run(session, "define axis fromPoint " + str(center.scene_coord[0]) + "," + str(center.scene_coord[1]) + "," + str(center.scene_coord[2]-100) + " toPoint " + str(center.scene_coord[0]) + "," + str(center.scene_coord[1]) + "," + str(center.scene_coord[2]+100))
        axis_X = run(session, "define axis fromPoint " + str(center.scene_coord[0]-100) + "," + str(center.scene_coord[1]) + "," + str(center.scene_coord[2]) + " toPoint " + str(center.scene_coord[0]+100) + "," + str(center.scene_coord[1]) + "," + str(center.scene_coord[2]))
        axis = run(session, f"define axis #{i+1:d}/A:348@ca #{i+1:d}/A:984@ca");
        shell_start = run(session, f"define centroid #{i+1:d}/A:348");
        bsol_center = run(session, f"define centroid #{i+1:d}/B:2940-3613");
        axis2 = run(session, "define axis fromPoint " + str(shell_start.scene_coord[0]) + "," + str(shell_start.scene_coord[1]) + "," + str(shell_start.scene_coord[2]) + " toPoint " + str(bsol_center.scene_coord[0]) + "," + str(bsol_center.scene_coord[1]) + "," + str(bsol_center.scene_coord[2]))
        flexion_tilt = (90-np.around(axis_Z[0].angle(axis[0]),2))*np.sign(axis[0].direction[2])*-1
        flexion_rot  = (90-np.around(axis_X[0].angle(axis[0]),2))*np.sign(axis[0].direction[0])
        shell2_tilt  = (90-np.around(axis_Z[0].angle(axis2[0]),2))*np.sign(axis2[0].direction[2])*-1

        axis_RY1     = run(session, f"define axis #{i+1:d}/A:915-932")
        axis_RY2     = run(session, f"define axis #{i+1:d}/A:979-1003")
        RY12_angle   = np.around(axis_RY1[0].angle(axis_RY2[0]),2)
        RY12d        = np.around(run(session, f"distance #{i+1:d}/A:915@ca #{i+1:d}/A:1003@ca"),2)

        distance1 = np.around(run(session, f"distance #{i+1:d}/A:4101@ca #{i+1:d}/D:4730@ca"),2)
        distance2 = np.around(run(session, f"distance #{i+1:d}/A:4075@ca #{i+1:d}/D:4736@ca"),2)

        pore_distance      = np.around(run(session, f"distance #{i+1:d}/A:4937@ca #{i+1:d}/C:4937@ca"),2)
        bsol_model_switch = run(session, "sel :3083")
        if len(bsol_model_switch.models)>0:
            RY12_BSol_distance = np.around(run(session, f"distance #{i+1:d}/A:979@ca #{i+1:d}/B:3083@ca"),2)
        else:
            RY12_BSol_distance = 0
        sel = run(session, "sel sequence TISPGDGRTF");cls=len(sel.models);
        f.write(files[i][0]+ "\t" + files[i][1] + f"\t{flexion_tilt:6.2f}\t{flexion_rot:6.2f}\t{shell2_tilt:6.2f}\t{cls:d}\t{RY12_angle:6.2f}\t{RY12d:6.2f}\t{distance1:6.2f}\t{distance2:6.2f}\t{pore_distance:6.2f}\t{RY12_BSol_distance:6.2f}\n")
        run(session, "close all")

