import glob
import os
import numpy as np
from chimerax.atomic import all_structures
from chimerax.core.commands import run

path = "models/"
path_out = '_situ_all'

files = glob.glob(path + "*.cif")

ali_name = ["TMD", "CSol", "S2S3", "Pore", "BSol"] # 5 rounds total
ali_ref = ['4540-4665,4787-4937',\
'3667-4059,4135-4174',\
'4666-4786',\
'4820-4956',\
'2939-3613']

ali_chain_id = [];
ali_chain_id.append(['/A,B,C,D:'])#TMD aligned on all chains
ali_chain_id.append(['/A:'])#CSol aligned on chain A 
ali_chain_id.append(['/D:'])#S2S3 aligned on chain D
ali_chain_id.append(['/A:'])#Pore aligned on chain A
ali_chain_id.append(['/B:'])#BSol aligned on chain B

target = [];target_name = [];target_chain_id = [];
#round 1; align on TMD and calculate shell1, shell2 and CSol on 1 chain. Here shell1 and shell2 should correlate with FlexionRot and FlexionTilt respectively, CSol should act as a pillow and correlate as well;
target_name.append(["CSol", "JBSol", "shell1"])
target.append(['3667-4059,4135-4174',\
'1657-2144,2145-2709,2955-3613',\
'1-1656'])
target_chain_id.append(['/A:'])

#round 2; align on CSol, TaF/CTD should correlate with Ca binding, S6c with pore opening, EF12 and S2S3 idk...
target_name.append(["EF12", "S2S3", "TaF", "S6c", "CTD"])
target.append(['4060-4134',\
'4666-4786',\
'4175-4253',\
'4938-4956',\
'4957-5037'])
target_chain_id.append(['/A:'])

#round 3, EF12 relative to the S2S3 of the neighbouring protomer, should correlate with states
target_name.append(["EF12"])
target.append(['4060-4134'])
target_chain_id.append(['/A:'])

#round 4, here to check if some membrane-mimicking systems actually have any significant impact, but don't expect much
target_name.append(["pVSD"])
target.append(['4540-4665,4787-4819'])
target_chain_id.append(['/A:'])

#round 5, here we check how RY1&2 moves relative to the post-RY3&4 BSol part...
target_name.append(["RY1&2"])
target.append(['850-1054'])
target_chain_id.append(['/A:'])

for a in range(len(ali_ref)):
    for t in range(len(target[a])):
        dir_name = "./scatters/" + target_name[a][t] + "_" + ali_name[a] + path_out
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        for i in range(len(files)):
            remain = files[i+1:]
            run(session, "open " + files[i])
            for j in range(len(remain)):
                if os.path.exists(dir_name + "/" + files[i][len(path):len(path)+4] + "_" + remain[j][len(path):len(path)+4] + ".txt"):
                    continue

                run(session, "open " + remain[j])
                run(session, "sel #1" + target_chain_id[a][0] + target[a][t] + " #2" + target_chain_id[a][0] +target[a][t])
                
                structures = all_structures(session)
                res0 = structures[0].residues; res1 = structures[1].residues
                num0 = set(res0[res0.selected].numbers); num1 = set(res1[res1.selected].numbers)
                nums = list(num1.intersection(num0))
                ranges = sum((list(t) for t in zip(nums, nums[1:]) if t[0]+1 != t[1]), [])
                iranges = iter(nums[0:1] + ranges + nums[-1:])
                sel = ','.join([str(n) + '-' + str(next(iranges)) for n in iranges])
              
                run(session, "align #2" + ali_chain_id[a][0] + ali_ref[a] + "@ca to #1" + ali_chain_id[a][0] + ali_ref[a] + "@ca matchNumbering 1 matchChainIds 1") # 1) align to reference
                align = run(session, "align #2"+ target_chain_id[a][0] + sel + "@ca toAtoms #1"+ target_chain_id[a][0] + sel + "@ca move nothing reportMatrix 1") # 2) align of target to fetch pre-aligned RMSD and rotation
                norm = np.round(align[2],6) #rmsd after domain-on-domain alignment for normalization
                rot = np.around(align[-1].axis_center_angle_shift()[-2],6) #domain rotation
                axis_shift = np.around(align[-1].axis_center_angle_shift()[-1],6) #shift along rotation axis
           
                rmsd = np.around(run(session, "rmsd #2" + target_chain_id[a][0] + sel + "@ca to #1" + target_chain_id[a][0] + sel + "@ca"),6) # RMSD of target, while aligned on reference
                shift = np.around(rmsd - norm,6) #rmsd-norm = displacement

                run(session, "close #2")
                with open(dir_name + "/" + files[i][len(path):len(path)+4] + "_" + remain[j][len(path):len(path)+4] + ".txt", 'w') as f:
    	             f.write(str(rot) + '\t' + str(axis_shift) + '\t' + str(rmsd) + '\t' + str(norm) + '\t' + str(shift) + '\t' + sel)
            run(session, "close #1")
