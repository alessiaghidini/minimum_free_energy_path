import numpy as np
import ast
import os
import subprocess
import configparser
import importlib

#################################################
####### check for all the modules install #######
required_modules = ["numpy", "ast", "os","subprocess", "configparser"]

missing_modules = []
for module in required_modules:
    try:
        importlib.import_module(module)
    except ImportError:
        missing_modules.append(module)

if missing_modules:
    print("Please install the following modules: {', '.join(missing_modules)}")
    print("You can install them using the following command:")
    print("pip install -r requirements.txt")
    exit(1)
#################################################
try:
    import MDAnalysis as mda
except ImportError:
    print("The MDAnalysis module is not installed. Try:")
    print("pip install MDAnalysis")
    exit(1)
#################################################

plumed_template_rmsd_file = 'plumed_rmsd_template.dat'
plumed_rmsd_file = 'plumed_rmsd.dat'
plumed_control_rmsd = 'colvar_rmsd'

mdp_template_file = "steered_md_template.mdp"
mdp_file = "steered_md_run.mdp"

plumed_template_file = 'plumed_template.dat'
plumed_run_file = 'plumed_run.dat'


####### read configuration file #######
def read_config():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_file_path = os.path.join(script_dir, 'input.ini')

    # Check if the file exists
    if not os.path.exists(config_file_path):
        print("Error: 'input.ini' file not found in the current directory.")
        exit()

    # Read the configuration file
    config = configparser.ConfigParser()
    try:
        config.read(config_file_path)
    except configparser.Error as e:
        print("Error reading 'input.ini' file: {}".format(e))
        exit()
   
    plumed = config.get("Paths", "plumed_path", fallback="")
    gromacs = config.get("Paths", "gromacs_path", fallback="")

    gromppCommand = config.get("Gromacs_Commands", "grompp_command", fallback="")
    mdrunCommand = config.get("Gromacs_Commands", "mdrun_command", fallback="")
    convertCommand = config.get("Gromacs_Commands", "convert_command", fallback="")

    number_frames = config.getint("Script_Parameters", "number_frames")
    equidistance_nm = config.getfloat("Script_Parameters", "equidistance_nm", fallback=1.0)
    tollerance_nm = config.getfloat("Script_Parameters", "tollerance_nm", fallback=1.0)

    KAPPA_rmsd1 = config.get("Plumed_Parameters", "KAPPA_rmsd1", fallback="")
    KAPPA_rmsd2 = config.get("Plumed_Parameters", "KAPPA_rmsd2", fallback="")


    selection_align = config.get("System_Info", "atoms_for_alignment", fallback="")
    selection_rmsd = config.get("System_Info", "atoms_for_rmsd", fallback="")
    entity0 = config.get("System_Info", "protein_ligand_atoms", fallback="")
    topology_file = config.get("System_Info", "topology_file", fallback="")
    index_file = config.get("System_Info", "index_file", fallback="")

    return (plumed,gromacs,
        gromppCommand, mdrunCommand, convertCommand,
        number_frames, equidistance_nm, tollerance_nm,
        KAPPA_rmsd1, KAPPA_rmsd2,
        selection_align,selection_rmsd,entity0,topology_file,index_file)


####### Plumed template rmsd file #######
def write_plumed_template_rmsd():
    with open(plumed_template_rmsd_file, 'w') as filout:
        filout.write('''rmsd1: RMSD REFERENCE=ref_pre TYPE=OPTIMAL
rmsd2: RMSD REFERENCE=ref_next TYPE=OPTIMAL
                     
PRINT ARG=rmsd1,rmsd2 STRIDE=1 FILE=colvar_rmsd

ENDPLUMED
Template to compute the RMSD between the structure and the previous one (rmsd1) and the next one (rmsd2).''')


####### Plumed run input file #######
def write_plumed_rmsd(ref_pre, ref_next):
    template = open(plumed_template_rmsd_file).read()
    out = open(plumed_rmsd_file, 'w')
    replacements = {'ref_pre':str(ref_pre), 'ref_next':str(ref_next)}
    for i in replacements.keys():
        template = template.replace(i, replacements[i])
    out.write(template)
    out.close()


####### calculate rmsd #######
def get_rmsd():
    colvar_rmsd = open(plumed_control_rmsd, 'r')
    line = colvar_rmsd.readlines()
    rmsd1 = float(line[1].split()[1])
    rmsd2 = float(line[1].split()[2])
    colvar_rmsd.close()
 
    return rmsd1, rmsd2


####### Plumed template input file #######
def write_plumed_template_input(entity0, KAPPA_rmsd1, KAPPA_rmsd2):
    with open(plumed_template_file, 'w') as filout:
        filout.write('''WHOLEMOLECULES ENTITY0={}
rmsd_1: RMSD REFERENCE=ref_pre TYPE=OPTIMAL
rmsd_2: RMSD REFERENCE=ref_next TYPE=OPTIMAL
MOVINGRESTRAINT ...
        ARG=rmsd_1,rmsd_2
        STEP0=0         AT0=rmsd1,rmsd2     KAPPA0=100.00,100.00
        STEP1=half_s    AT1={},rmsd2    KAPPA1={},{}
        STEP2=step      AT2={},rmsd2    KAPPA2={},{}
        STEP3=step_t    AT3={},rmsd2    KAPPA3={},{}
... MOVINGRESTRAINT
PRINT ARG=rmsd_1,rmsd_2 FILE=colvar STRIDE=10 FMT=%8.4f'''.format(entity0, equidistance_nm, KAPPA_rmsd1,KAPPA_rmsd2, (equidistance_nm-tollerance_nm) ,KAPPA_rmsd1,float(KAPPA_rmsd2)/2, (equidistance_nm-tollerance_nm) ,KAPPA_rmsd1,float(KAPPA_rmsd2)/5))

####### Plumed run input file #######
def write_plumed_input(steps_md, ref_pre, ref_next, rmsd_1, rmsd_2):
    template = open(plumed_template_file).read()
    out = open(plumed_run_file, 'w')
    replacements = {'half_s':str(steps_md/2), 'step':str(steps_md), 'step_t':str(steps_md+10000), 'ref_pre':str(ref_pre), 'ref_next':str(ref_next), 'rmsd1':str(rmsd_1), 'rmsd2':str(rmsd_2)}
    for i in replacements.keys():
        template = template.replace(i, replacements[i])
    out.write(template)
    out.close()


####### mdp gromacs input file #######
def write_gromacs_input(steps_md):
    template = open(mdp_template_file).read()
    out = open(mdp_file, 'w')
    replacements = {'numsteps':str(steps_md)}
    for i in replacements.keys():
        template = template.replace(i, replacements[i])
    out.write(template)
    out.close()


########## define nearest structure ########## 
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))    


########## select alignment and rmsd ######## 
def write_reference_pdb(gro_file, output_pdb):
    u = mda.Universe(gro_file)

    align = u.select_atoms(selection_align)
    rmsd  = u.select_atoms(selection_rmsd)

    u.add_TopologyAttr("occupancies")
    u.add_TopologyAttr("tempfactors")

    u.atoms.occupancies = 0.0
    u.atoms.tempfactors = 0.0

    u.atoms[align.indices].occupancies = 1.0
    u.atoms[rmsd.indices].tempfactors = 1.0

    u.atoms.write(output_pdb)

    return output_pdb

def filter_pdb_by_flags(input_pdb, output_pdb):
    with open(input_pdb, "r") as f:
        lines = f.readlines()

    kept = []

    for line in lines:
        if not line.startswith("ATOM"):
            kept.append(line)
            continue

        # PDB fixed columns:
        occupancy = float(line[54:60])
        beta = float(line[60:66])

        # keep if selected in either mask
        if occupancy == 1.0 or beta == 1.0:
            kept.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(kept)

    return output_pdb

##########################################

# read input
(plumed, gromacs, gromppCommand, mdrunCommand, convertCommand,
 number_frames, equidistance_nm, tollerance_nm, KAPPA_rmsd1, KAPPA_rmsd2,
 selection_align, selection_rmsd, entity0, topology_file, index_file) = read_config()


# prepare log file
log_file = open('equidistant_frames.log', 'w')
log_file.write('PROGRAM TO MAKE EQUIDISTANT THE FRAMES STARTING\n\n')

log_file_index_list = open('gro_for_path.log', 'w')

# prepare template plumed input to compute rmsd and perform steered md
write_plumed_template_rmsd()
write_plumed_template_input(entity0, KAPPA_rmsd1,KAPPA_rmsd2)

# read input file
#plumed_path,gromacs_path,topology_file,index_file,equidistance_nm,KAPPA_rmsd1,KAPPA_rmsd2 = read_config(config_file)

# final list of pdb for final path
pdb_final_path = []

# str_steered: structure to be moved during steered md
# str_pre: previous structure (ref to compute rmsd_1)
# str_next: following structure (ref to compute rmsd_2) 

list_frames = list(range(1,(number_frames + 1)))
frame_pre = list_frames[0]
frame_steered = frame_pre

pdb_final_path.append(f"{frame_pre}.pdb")
log_file_index_list.write(f"{frame_pre}.pdb, ")

# create pdb of first frame for path 
write_reference_pdb(f"{frame_pre}.gro",f"{frame_pre}.pdb")
filter_pdb_by_flags(f"{frame_pre}.pdb",f"{frame_pre}.pdb")


# equispacing is repeated until frame_steered of previous replica is lower than list_frames-1
# i.e., until new frame to steered is lower than list_frames
while frame_steered < len(list_frames)-1:    
    # define structure to move and refs
    frame_steered = frame_steered+1
    frame_next = frame_steered+1
    
    gro_steered = f"{frame_steered}.gro"
    ref_pre = f"{frame_pre}_ref.pdb"
    ref_next = f"{frame_next}_ref.pdb"

    # create reference pdb
    # previous
    write_reference_pdb(f"{frame_pre}.gro",ref_pre)
    filter_pdb_by_flags(ref_pre,ref_pre)
    # next
    write_reference_pdb(f"{frame_next}.gro",ref_next)
    filter_pdb_by_flags(ref_next,ref_next)
    
    # compute rmsd1 with ref_pre and rmsd2 with ref_next
    write_plumed_rmsd(ref_pre, ref_next)
    plumed_rmsd = subprocess.Popen(plumed + "plumed driver --plumed %s --igro %s "%(plumed_rmsd_file,gro_steered), shell=True)
    plumed_rmsd.wait()
    rmsd1, rmsd2 = get_rmsd()
    log_file.write('\n\nstarting iteration to move %d'%(frame_steered))
    log_file.write('\nrmsd1:{%f} - between {%s} and previous reference {%s} \nrmsd2:{%f} - between {%s} and next reference {%s}'%(rmsd1,frame_steered,ref_pre,rmsd2,frame_steered,ref_next))
    log_file.flush()

    while (rmsd1 + tollerance_nm) < equidistance_nm:
        frame_steered = frame_next
        frame_next = frame_steered+1
        gro_steered = f"{frame_steered}.gro"
        ref_next = f"{frame_next}_ref.pdb"
        # next
        write_reference_pdb(f"{frame_next}.gro",ref_next)
        filter_pdb_by_flags(ref_next,ref_next)

        # compute rmsd1 with ref_pre and rmsd2 with ref_next
        write_plumed_rmsd(ref_pre, ref_next)
        plumed_rmsd = subprocess.Popen(plumed + "plumed driver --plumed %s --igro %s "%(plumed_rmsd_file,gro_steered), shell=True)
        plumed_rmsd.wait()
        rmsd1, rmsd2 = get_rmsd()
        log_file.write('\nrmsd1<=equidistance_nm')
        log_file.write('\n\nstarting iteration to move %d'%(frame_steered))
        log_file.write('\nrmsd1:{%f} - between {%d} and previous reference {%s} \nrmsd2:{%f} - between {%d} and next reference {%s}'%(rmsd1,frame_steered,ref_pre,rmsd2,frame_steered,ref_next))
        log_file.flush()

    
    # estimee of intermediate frames, thus md steps for steered md 
    interm_frames = int(rmsd1/equidistance_nm)
    log_file.write('\nexpected number of intermediate frames %d'%(interm_frames))

    # steered md steps to perform (100 ps for intermediate frame) and inputs
    steps_md = interm_frames*50000
    write_gromacs_input(steps_md)
    write_plumed_input(steps_md, ref_pre, ref_next, rmsd1, rmsd2)

    while interm_frames>= 1:
        
        # perform steered md for all intermediate frames
        log_file.write('\nsteered md for {%d} - {%d} intermediate frames to be moved'%(frame_steered,interm_frames))
        grompp_process = subprocess.Popen('{}{} -f {} -c {} -r {} -p {} -n {} -maxwarn 5 -o steered_md.tpr'.format(gromacs, gromppCommand, mdp_file, gro_steered, gro_steered, topology_file, index_file),shell=True)
        grompp_process.wait()
        mdrun_process = subprocess.Popen('{}{} -deffnm steered_md -plumed {} -ntomp 8'.format(gromacs, mdrunCommand, plumed_run_file),shell=True)
        mdrun_process.wait()
        log_file.write('\n finished steered md for ' + str(frame_steered))

        # check to see if there is a frame satisfing the equidistance condition
        plumed_rmsd = subprocess.Popen(plumed + "plumed driver --plumed "+plumed_rmsd_file+" --mf_xtc steered_md.xtc", shell=True)
        plumed_rmsd.wait()
        f=open(plumed_control_rmsd,"r")
        lines=f.readlines()[1:]
        rmsd1_list=[]
        for j in lines:
            rmsd1_list.append(float(j.split()[1]))
        f.close()
        rmsdmin = nearest(rmsd1_list, equidistance_nm)

        # if condition is satisfied by a frame:
        if (np.abs(rmsdmin - equidistance_nm) < tollerance_nm):
            
            log_file.write('\n structure %s_%d after steered md satisfies equidistant condition'%(frame_steered,interm_frames))
            log_file.write('\n adding structure %s_%d to path'%(frame_steered,interm_frames))
            frame_equidistant = rmsd1_list.index(nearest(rmsd1_list, equidistance_nm))
            
            # load frame with correct equidistance and save to gro file
            new_str = mda.Universe(f"{frame_steered}.gro", "steered_md.xtc")
            new_str.trajectory[frame_equidistant]
            new_str.atoms.write(f"{frame_steered}_{interm_frames}.gro")
            
            # frame save as pdb frame and add to path  
            # this becomes the new previous reference (next reference is the same and structure steered is the same)
            frame_next = f"{frame_steered}_{interm_frames}"
            ref_next = f"{frame_next}_ref.pdb"
            pdb_final_path.append(ref_next)
            log_file_index_list.write(f"{frame_next}_ref.pdb, ")
            # next
            write_reference_pdb(f"{frame_next}.gro",ref_next)
            filter_pdb_by_flags(ref_next,ref_next)
            
            # compute rmsd1 with ref_pre and rmsd2 with ref_next
            write_plumed_rmsd(ref_pre, ref_next)
            rmsd1, rmsd2 = get_rmsd()

            #update of new intermediate frames
            interm_frames -= 1

            # steered md steps to perform (100 ps for intermediate frame) and inputs
            steps_md = interm_frames*50000
            write_gromacs_input(steps_md)
            write_plumed_input(steps_md, ref_pre, ref_next, rmsd1, rmsd2)

            subprocess.call("rm bc* \#*", shell=True)
            log_file.flush()
        
        # if condition is NOT satisfied by a frame:
        else:

            # the simulation is extended until the condition is satisfied          
                while (True):
                    log_file.write('\nextension of steered md for %s - %d intermediate frames to be moved'%(frame_steered,interm_frames))
                    convert_process = subprocess.Popen('{}{} -s steered_md.tpr -extend 50 -o steered_md.tpr'.format(gromacs,convertCommand), shell=True)
                    convert_process.wait()
    
                    mdrun_process = subprocess.Popen('{}{} -deffnm steered_md -cpi steered_md.cpt -plumed {} -ntomp 8'.format(gromacs, mdrunCommand, plumed_run_file),shell=True)
                    mdrun_process.wait()
    
                    log_file.write('\n finished steered md for ' + str(frame_steered))

                    # check to see if there is a frame satisfing the equidistance condition
                    plumed_rmsd = subprocess.Popen(plumed + "plumed driver --plumed %s --mf_xtc steered_md.xtc"%(plumed_rmsd_file), shell=True)
                    plumed_rmsd.wait()
                    f=open(plumed_control_rmsd,"r")
                    lines=f.readlines()[1:]
                    rmsd1_list=[]
                    for j in lines:
                        rmsd1_list.append(float(j.split()[1]))
                    f.close()
                    rmsdmin = nearest(rmsd1_list, equidistance_nm)
                    
                    # if condition is satisfied by a frame:
                    if (np.abs(rmsdmin - equidistance_nm) <= tollerance_nm):
                        # the frame is added to the path
                        log_file.write('\n structure %s_%d after steered md satisfies equidistant condition'%(frame_steered,interm_frames))
                        log_file.write('\n adding structure %s_%d to path'%(frame_steered,interm_frames))
                        frame_equidistant = rmsd1_list.index(nearest(rmsd1_list, equidistance_nm))


                        # load frame with correct equidistance and save to gro file
                        new_str = mda.Universe(f"{frame_steered}.gro", "steered_md.xtc")
                        new_str.trajectory[frame_equidistant]
                        new_str.atoms.write(f"{frame_steered}_{interm_frames}.gro")
            
                        # frame save as pdb frame and add to path  
                        # this becomes the new previous reference (next reference is the same and structure steered is the same)
                        frame_next = f"{frame_steered}_{interm_frames}"
                        ref_next = f"{frame_next}_ref.pdb"
                        pdb_final_path.append(ref_next)
                        log_file_index_list.write(f"{frame_next}_ref.pdb, ")
                        # next
                        write_reference_pdb(f"{frame_next}.gro",ref_next)
                        filter_pdb_by_flags(ref_next,ref_next)

                        # compute rmsd1 with ref_pre and rmsd2 with ref_next
                        write_plumed_rmsd(ref_pre, ref_next)
                        rmsd1, rmsd2 = get_rmsd()

                        #update of new intermediate frames
                        interm_frames -= 1
                        #log_file.write('\n\nstarting iteration %d_%d'%(frame_steered,interm_frames))

                        # steered md steps to perform (100 ps for intermediate frame) and inputs
                        steps_md = interm_frames*50000
                        write_gromacs_input(steps_md)
                        write_plumed_input(steps_md, ref_pre, ref_next, rmsd1, rmsd2)

                        # restart with new previous reference
                        subprocess.call("rm bc* \#*", shell=True)
                        #log_file.write('\nrmsd1:{%f} - between {%s} and previous reference {%s} \nrmsd2:{%f} - between {%s} and next reference {%s}'%(rmsd1,frame_steered,ref_pre,rmsd2,frame_steered,ref_next))
                        break

    interm_frames = 1
    frame_pre = str(frame_steered)+"_"+str(interm_frames)

# save the final path as pdb
    
first = mda.Universe(pdb_final_path[0])
with mda.Writer("reference_path.pdb", first.atoms.n_atoms) as W:
    for pdb in pdb_final_path:
        u = mda.Universe(pdb)
        W.write(u.atoms)    

log_file.write('\n\n Program finished - path_reference.pdb created')
log_file.flush()
