#!/usr/bin/env python3

import os,sys
import numpy as np
import time
import multiprocessing
import re
import subprocess
import logging

#########
#
# Description :
# 	run all namelist.ssh files and all pythonscript for graphs
#
# Last update : 30/09/2024
#
# Authors: Zhizhao Wang and Youngseob Kim
#########

# == User configuration
arg_run = True
arg_figure = True
multiproc_run = True

# Choice of cases
# 0: all
# 1: condensation
# 2: coagulation
# 3: nucleation
# 4: viscosity
# 5: caco3
# 6: mono-terpene
# 7: mcm
# 8: user-defined 
case_type = 0

# if use_current_directory = True,
# input_directory and output_directory are ignored. 
# use_current_directory = True
# input_directory = "/cerea_raid/users/kimy/work/ssh-aerosol/ssh-aerosol.git/"
# output_directory = "/net/libre/merida/kimy/ssh-aerosol-testcase/"
# ==

if case_type == 1:
    keyword = 'cond'
elif case_type == 2:
    keyword = 'coag'
elif case_type == 3:
    keyword = 'nucl'
elif case_type == 4:
    keyword = 'visc'
elif case_type == 5:  
    keyword = 'caco3'
elif case_type == 6:
    keyword = 'prt4act3'
elif case_type == 7:  
    keyword = 'monoterpene'
elif case_type == 8:  
    keyword = 'mcm'
elif case_type == 9:  
    keyword = 'vocox'    

"""
Get date and time as a string
"""
def get_now():        

    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    return dt_string


def write_status(status, cmd):

    if (status == 0):
        message = ' => Success'
    else:
        message = ' => Fail'
    logging.info(cmd.ljust(50) + message.ljust(20))


def run_monoterpene_cases():

    ssh_dir = "./"
    cwd = os.getcwd()


    # SOA formation from monoterpene oxidation

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -g=fast"
    logging.getLogger(cmd)
    status = os.system(cmd)
    write_status(status, cmd)
    
    namelist = "namelist_mt_ref.ssh"
    print("Run: namelist file : " + namelist)
    logging.getLogger("Run: namelist file : " + namelist)
    
    cmd = "./ssh-aerosol INIT/" + namelist + " 1" + \
        " > results/ssh_log_" + namelist + " 2>&1"
    status = os.system(cmd)
    write_status(status, cmd)
    
    namelist = "namelist_mt_rdc.ssh"
    print("Run: namelist file : " + namelist)
    logging.getLogger("Run: namelist file : " + namelist)
    
    cmd = "./ssh-aerosol INIT/" + namelist + " 1" + \
        " > results/ssh_log_" + namelist + " 2>&1"
    status = os.system(cmd)
    write_status(status, cmd)    

    os.chdir(ssh_dir)

def run_mcm_cases():
    
    os.system("bash INIT/launch_mcm_diffRO2.sh")
    os.system("bash INIT/launch_mcm_diffRO2_smiles.sh")

def run_case(namelist):

    cmd = "./ssh-aerosol INIT/" + namelist + \
        " > results/ssh_log_" + namelist + " 2>&1"
    status = os.system(cmd)
    write_status(status, cmd)   
    
def run_user_defined_cases():

    ssh_dir = "./"
    cwd = os.getcwd()

    # Run the base case

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme.cfg -m=species_matching.dat"
    logging.getLogger(cmd)
    status = os.system(cmd)
    write_status(status, cmd)

    namelist = "namelist_vocox.ssh"
    run_case(namelist)

    namelist = "namelist_vocox_dNO2.ssh"
    run_case(namelist)    

    namelist = "namelist_vocox_mNO2.ssh"
    run_case(namelist)
    
    # Run the EXPL case

    print ("=== Run EXPL case ===")
    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_expl.cfg -m=species_matching_expl.dat"
    logging.getLogger(cmd)
    status = os.system(cmd)
    write_status(status, cmd)   

    namelist = "namelist_vocox_expl.ssh"
    run_case(namelist)

    namelist = "namelist_vocox_expl_dNO2.ssh"
    run_case(namelist)

    namelist = "namelist_vocox_expl_mNO2.ssh"
    run_case(namelist)

    # Run the RDC case

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_rdc.cfg -m=species_matching_rdc.dat"
    logging.getLogger(cmd)
    status = os.system(cmd)
    write_status(status, cmd)   

    namelist = "namelist_vocox_rdc.ssh"
    run_case(namelist)

    namelist = "namelist_vocox_rdc_dNO2.ssh"
    run_case(namelist)

    namelist = "namelist_vocox_rdc_mNO2.ssh"
    run_case(namelist)
    
    # Run the noautox case

    print ("=== Run noautox case ===")
    
    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_h2o-noautox.cfg -m=species_matching_h2o-noautox.dat"
    logging.getLogger(cmd)
    status = os.system(cmd)
    write_status(status, cmd)   

    namelist = "namelist_vocox-h2onoautox.ssh"
    run_case(namelist)

    namelist = "namelist_vocox-h2onoautox_dNO2.ssh"
    run_case(namelist)

    namelist = "namelist_vocox-h2onoautox_mNO2.ssh"
    run_case(namelist)
    
    os.chdir(cwd)


def multiprocessing_func(k):

    status = -999
    
    if k[-4:] == '.ssh':

        print("Run: namelist file : " + k)
        logging.getLogger("Run: namelist file : " + k)
            
        # run simulations
        cmd = "time ./ssh-aerosol INIT/" + k + \
            " > results/ssh_log_" + k[:-4] + " 2>&1"
        status = os.system(cmd)

    if (status != 0):
        sys.exit(1)
        
    return status

def run_prt4act3_case():
    """This function runs the test case with namelist_prt4act3.ssh"""

    # Check if files exist
    filename = "INIT/namelist_prt4act3.ssh"
    shellname = "INIT/launch_prt4act3.sh"
    for f in [filename, shellname]:
        if not os.path.exists(f):
            logging.error(f"{f} does not exist.")
            return

    # Run the test case
    print(f"\nRunning {shellname} ...")

    cmd = f"time ./{shellname}"
    status = os.system(cmd)
    if status != 0:
        logging.error(f"Failed to run the test case with {shellname}")
        return


if __name__ == '__main__':

    starttime = time.time()
    processes = []
    processes_name = []
    simulation_dir = './'
    
    # Make logfile
    logfile = "run_testcase.log." + get_now() 
    try:
        os.remove(logfile)
        print("Remove the existing logfile " + logfile)
    except OSError:
        pass
    
    logging.basicConfig(filename = simulation_dir + "/" + logfile,
                        level = logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    console = logging.StreamHandler()
    console.setLevel(logging.ERROR)
    logging.getLogger("").addHandler(console)

        
    ssh_dir = simulation_dir
    os.chdir(ssh_dir)
       
    if (arg_run):

        cmd = "clean -a"
        os.system(cmd)

        cmd = "compile"
        logging.getLogger(cmd)
        status = os.system(cmd)
        write_status(status, cmd)   
        
        # remove all result files.
        subprocess.run(['rm', '-rf', 'results/*'])

        # get all file names in the folder INIT
        if case_type == 0:
            initfile = os.listdir(ssh_dir+ 'INIT/')
        else:
            initfile = [f for f in os.listdir(ssh_dir+ 'INIT/') if re.search(keyword,f)]

        for i in range(len(initfile)):
            if (initfile[i][-4:]==".ssh"):

                case_name = initfile[i] 

                # Skip some test cases. They are launched after.
                if (case_name == "namelist_mt_ref.ssh" or \
                    case_name == "namelist_mt_rdc.ssh" or \
                    case_name == "namelist_mcm-wFGL.ssh" or \
                    case_name == "namelist_prt4act3.ssh" or \
                    case_name == "namelist_mcm-wSMILES.ssh" or \
                    case_name[9:14] == "vocox"):
                    continue
                    
                logging.info("run the case " + initfile[i])
                        
                if (multiproc_run):
                    p = multiprocessing.Process(target = multiprocessing_func, args = (initfile[i],))
                    processes.append(p)
                    processes_name.append(initfile[i])
                    p.start()
                else:
                    status = multiprocessing_func(initfile[i])
                    if (status == 0):
                        message = ' => Success'
                    else:
                        message = ' => Fail'
                    logging.info(initfile[i].ljust(50) + message.ljust(20))

        if (multiproc_run):                    
            for process, process_name in zip(processes, processes_name):
                process.join()
                status = process.exitcode
                if (status == 0):
                    message = ' => Success'
                else:
                    message = ' => Fail'
                logging.info(process_name.ljust(50) + message.ljust(20))
        
        # Rub prt4act3 case
        if case_type == 0 or case_type == 6:
            run_prt4act3_case()
        
        # Run monoterpene cases
        if case_type == 0 or case_type == 7:
            run_monoterpene_cases()

        # Run MCM cases
        if case_type == 0 or case_type == 8:                
            run_mcm_cases()

        # Run user-defined cases
        if case_type == 0 or case_type == 9:
            run_user_defined_cases()

        logging.info('That took {} seconds'.format(time.time() - starttime))

        logging.info('! ! All simulations are done, ')

        logging.info('! Run all postprocessing python scripts in graph folder: ')

    if (arg_figure):
        # record the number of python scripts
        num = 0

        if case_type == 0:
            # get all file names in the folder graph
            graphfile =  os.listdir(ssh_dir+ 'graph/')
        else:
            graphfile = [f for f in os.listdir(ssh_dir+ 'graph/') if re.search(keyword,f)]                

        os.chdir('graph')
        for k in graphfile :

            if k[-3:] == '.py':

                num = num + 1
                logging.info("Run "+str(num) + ", python script : " + k)
                # run python script
                cmd = 'python3 ' + k
                status = os.system(cmd)
                write_status(status, cmd)   
                
        logging.info('! ! All scripts are done, ' + 'number of .py files : ' + str(num))


            
