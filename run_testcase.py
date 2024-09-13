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
# Last update : 15/03/2022
#
# Authors: Zhizhao Wang and Youngseob Kim
#########

# == User configuration
arg_run = True
arg_figure = True
multiproc_run = True

# if use_current_directory = True,
# input_directory and output_directory are ignored. 
use_current_directory = True
input_directory = "/cerea_raid/users/kimy/work/ssh-aerosol/ssh-aerosol.git/"
output_directory = "/net/libre/merida/kimy/ssh-aerosol-testcase/"
# ==

"""
Get date and time as a string
"""
def get_now():        

    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    return dt_string


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

        if (use_current_directory):
            simulation_dir = './'
        else:
            simulation_time = get_now()

            simulation_dir = output_directory + "/" + simulation_time + "/"
            os.makedirs(simulation_dir + "/results")
            os.makedirs(simulation_dir + "/graph")

            # Copy files
            infile_list = ["ssh-aerosol",
                           "src",
                           "INIT",
                           "coef*",
                           "graph",
                           "inputs",
                           "photolysis",
                           "species-list"]
        
            for infile in infile_list:
                os.system("cp -Lrf " + input_directory + "/" + 
                          infile + " " + simulation_dir + "/")
            
            print("Make result directory at ", simulation_dir)
            

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
