import subprocess
import os
import time
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def call(cmd, shell=False):
    start = time.time()
    if shell:
        print(">>> {}".format(cmd))
        os.system(cmd)
    else:
        print(">>> {}".format(cmd))
        cmd_list = re.split("\\s+", cmd)
        p = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            print(bcolors.FAIL + err + bcolors.FAIL)
        if out:
            print(bcolors.OKGREEN + out + bcolors.OKGREEN)

        if p.returncode != 0:
            print(bcolors.FAIL + "EXTERNAL CALL FAILED. TERMINATING..." + bcolors.FAIL)
            exit(1)

    end = time.time()
    print(">>> {}. FINISHED in {:.0f}s".format(cmd, end - start))
