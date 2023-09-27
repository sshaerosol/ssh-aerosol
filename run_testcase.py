import os,sys
import numpy as np
import time
import multiprocessing
import re
import subprocess

#########
#
# Description :
# 	run all namelist.ssh files and all pythonscript for graphs
#
# Last update : 15/03/2022
#
# Authors: Zhizhao Wang and Youngseob Kim
#########

arg_run = True
arg_figure = True
multiproc_run = True

def multiprocessing_func(k):
        if k[-4:] == '.ssh':
            print("Run: namelist file : " + k)
            # run simulations
            cmd = "time ./ssh-aerosol INIT/" + k + \
            " > results/ssh_log_" + k[:-4] + " 2>&1"
            os.system(cmd)

if __name__ == '__main__':
        starttime = time.time()
        processes = []

        ssh_dir = './'
        os.chdir(ssh_dir)

        if (os.path.isfile('./ssh-aerosol') == False):
            print("IOError: did you compile the program?")
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
                        print("run the case ", initfile[i])

                        if (multiproc_run):
                                p = multiprocessing.Process(target = multiprocessing_func, args = (initfile[i],))
                                processes.append(p)
                                p.start()
                        else:
                                multiprocessing_func(initfile[i])

            if (multiproc_run):                    
                for process in processes:
                    process.join()

            print('That took {} seconds'.format(time.time() - starttime))

            print('! ! All simulations are done, ')

            print('! Run all postprocessing python scripts in graph folder: ')

        if (arg_figure):
            # record the number of python scripts
            num = 0
            # get all file names in the folder graph
            graphfile =  os.listdir(ssh_dir+ 'graph/')
            os.chdir('graph')
            for k in graphfile :
	            if k[-3:] == '.py':
	                num = num + 1
	                print("Run "+str(num) + ", python script : " + k)
                        # run python script
	                os.system('python3 ' + k)
            print('! ! All processing file are done, ' + 'number of .py files : ' + str(num))
