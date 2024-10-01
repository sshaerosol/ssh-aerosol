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

# if use_current_directory = True,
# input_directory and output_directory are ignored. 
# use_current_directory = True
# input_directory = "/cerea_raid/users/kimy/work/ssh-aerosol/ssh-aerosol.git/"
# output_directory = "/net/libre/merida/kimy/ssh-aerosol-testcase/"
# ==

"""
Get date and time as a string
"""
def get_now():        

    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    return dt_string


def run_monoterpene_cases():

    ssh_dir = "./"
    cwd = os.getcwd()


    # SOA formation from monoterpene oxidation

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -g=fast"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_mt_ref.ssh 1"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_mt_rdc.ssh 1"
    os.system(cmd)

    os.chdir(ssh_dir)

def run_mcm_cases():
    
    os.system("bash INIT/launch_mcm_diffRO2.sh")
    os.system("bash INIT/launch_mcm_diffRO2_smiles.sh")

def run_user_defined_cases():

    ssh_dir = "./"
    cwd = os.getcwd()

    # Run the base case

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme.cfg -m=species_matching.dat"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_dNO2.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_mNO2.ssh"
    os.system(cmd)

    # Run the EXPL case

    print ("=== Run EXPL case ===")
    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_expl.cfg -m=species_matching_expl.dat"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_expl.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_expl_dNO2.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_expl_mNO2.ssh"
    os.system(cmd)


    # Run the RDC case

    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_rdc.cfg -m=species_matching_rdc.dat"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_rdc.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_rdc_dNO2.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox_rdc_mNO2.ssh"
    os.system(cmd)

    # Run the noautox case

    print ("=== Run noautox case ===")
    
    cmd = "clean -a"
    os.system(cmd)

    cmd = "compile -b=user -c=cb05-ozone -u=user_defined_scheme_h2o-noautox.cfg -m=species_matching_h2o-autox.dat"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox-h2onoautox.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox-h2onoautox_dNO2.ssh"
    os.system(cmd)

    cmd = "./ssh-aerosol INIT/namelist_vocox-h2onoautox_mNO2.ssh"
    os.system(cmd)
    
    os.chdir(cwd)


def multiprocessing_func(k):
        if k[-4:] == '.ssh':
            print("Run: namelist file : " + k)
            logging.getLogger("Run: namelist file : " + k)
            
            # run simulations
            cmd = "time ./ssh-aerosol INIT/" + k + \
            " > results/ssh_log_" + k[:-4] + " 2>&1"
            os.system(cmd)

if __name__ == '__main__':
        starttime = time.time()
        processes = []

        simulation_dir = './'

        # Make logfile    
        logging.basicConfig(filename = simulation_dir + '/run_testcase.log',
                            level = logging.DEBUG,
                            format='%(asctime)s %(message)s',
                            datefmt='%d/%m/%Y %H:%M:%S')
        console = logging.StreamHandler()
        console.setLevel(logging.ERROR)
        logging.getLogger("").addHandler(console)

        
        ssh_dir = simulation_dir
        os.chdir(ssh_dir)

        if (os.path.isfile('./ssh-aerosol') == False):
            logging.info("IOError: did you compile the program?")
            sys.exit()
        
        if (arg_run):
            # remove all result files.
            subprocess.run(['rm', '-rf', 'results/*'])

            # get all file names in the folder INIT
            all_files = True
            if all_files:
                initfile = os.listdir(ssh_dir+ 'INIT/')
            else:
                keyword = 'monomer' # 'soalp'
                initfile = [f for f in os.listdir(ssh_dir+ 'INIT/') if re.search(keyword,f)]

            for i in range(len(initfile)):
                if (initfile[i][-4:]==".ssh"):

                    case_name = initfile[i] 

                    # Skip some test cases. They are launched after.
                    if (case_name == "namelist_mt_ref.ssh" or \
                        case_name == "namelist_mt_rdc.ssh" or \
                        case_name == "namelist_mcm-wFGL.ssh" or \
                        case_name == "namelist_mcm-wSMILES.ssh" or \
                        case_name[9:14] == "vocox"):
                        continue
                    
                    logging.info("run the case " + initfile[i])
                        
                    if (multiproc_run):
                        p = multiprocessing.Process(target = multiprocessing_func, args = (initfile[i],))
                        processes.append(p)
                        p.start()
                    else:
                        multiprocessing_func(initfile[i])

            if (multiproc_run):                    
                for process in processes:
                    process.join()

            # Run monoterpene cases
            run_monoterpene_cases()

            # Run MCM cases
            run_mcm_cases()

            # Run user-defined cases
            run_user_defined_cases()

            logging.info('That took {} seconds'.format(time.time() - starttime))

            logging.info('! ! All simulations are done, ')

            logging.info('! Run all postprocessing python scripts in graph folder: ')

        if (arg_figure):
            # record the number of python scripts
            num = 0
            # get all file names in the folder graph
            graphfile =  os.listdir(ssh_dir+ 'graph/')
            os.chdir('graph')
            for k in graphfile :
	            if k[-3:] == '.py':

	                num = num + 1
	                logging.info("Run "+str(num) + ", python script : " + k)
                        # run python script
	                os.system('python3 ' + k)
            logging.info('! ! All processing file are done, ' + 'number of .py files : ' + str(num))


            
