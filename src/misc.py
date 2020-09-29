import subprocess
import os
import time
import re

class bcolors:
    FAIL = '\033[47m\033[91m'
    ERROR = '\033[96m'
    COMMAND = '\033[4m'
    WARNING = '\033[33m'
    OUTPUT = ''

def call(cmd, shell=False):
    start = time.time()
    if shell:
        print(">>> {}".format(cmd))
        os.system(cmd)
        end = time.time()
        print(">>> {}. FINISHED in {:.0f}s".format(cmd, end - start))
    else:
        print(">>> {1}{0}{1}".format(cmd, bcolors.COMMAND))
        cmd_list = re.split("\\s+", cmd)
        p = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.getcwd())
        out, err = p.communicate()

        if err:
            err_list = str(err).splitlines()
        else:
            err_list = []

        if out:
            out_list = str(out).splitlines()
        else:
            out_list = []

        end = time.time()
        print("FINISHED in {:.0f}s".format(end - start))

        if err:
            for l in err_list:
                print("{1}{0}{1}".format(l, bcolors.WARNING))
        if out:
            for l in out_list:
                print("{1}{0}{1}".format(l, bcolors.OUTPUT))

        if p.returncode != 0:
            print(bcolors.FAIL + "EXTERNAL CALL FAILED. TERMINATING..." + bcolors.FAIL)
            exit(1)
