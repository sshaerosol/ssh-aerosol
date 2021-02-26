import os,sys
import numpy as np
import time
import multiprocessing

#########
#
# Description :
# 	run all namelist.ssh files and all pythonscript for graphs
#
# Last update : 05/07/2019
#
# Authors: Zhizhao Wang and Youngseob Kim
#########

def multiprocessing_func(k):
        if k[-4:] == '.ssh':
	        print("Run: namelist file : " + k)
                # run simulations
	        os.system('time ssh-aerosol INIT/' + k )

if __name__ == '__main__':
        starttime = time.time()
        processes = []

        ssh_dir = './'
        os.chdir(ssh_dir)

        # remove all result files.
        os.system('rm -rf results')

        # get all file names in the folder INIT
        initfile = os.listdir(ssh_dir+ 'INIT/')
        for i in range(len(initfile)):
                p = multiprocessing.Process(target = multiprocessing_func, args = (initfile[i],))
                processes.append(p)
                p.start()

        for process in processes:
                process.join()

        print('That took {} seconds'.format(time.time() - starttime))

        print('! ! All simulations are done, ')

        print('! Run all postprocessing python scripts in graph folder: ')

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
	                os.system('python ' + k)
        print('! ! All processing file are done, ' + 'number of .py files : ' + str(num))
