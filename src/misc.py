import subprocess
import os
import time
import re

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class bcolors:
    FAIL = '\033[47m\033[91m'
    ERROR = '\033[96m'
    COMMAND = '\033[4m'
    WARNING = '\033[33m'
    OUTPUT = ''
    END = '\033[0m'

def call(cmd, shell=False):
    start = time.time()
    if shell:
        print(">>> {}{}{}".format(bcolors.COMMAND, cmd, bcolors.END))
        r = os.system(cmd)
        end = time.time()

        if r != 0:
            print("{}EXTERNAL CALL FAILED (after {:.0f}s). TERMINATING...{}".format(bcolors.FAIL, end - start, bcolors.END))
            exit(1)

        print("FINISHED in {:.0f}s".format(end - start))
    else:
        print(">>> {}{}{}".format(bcolors.COMMAND, cmd, bcolors.END))
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
                print("{}{}{}".format(bcolors.WARNING, l, bcolors.END))
        if out:
            for l in out_list:
                print("{}{}{}".format(bcolors.OUTPUT, l, bcolors.END))

        if p.returncode != 0:
            print(bcolors.FAIL + "EXTERNAL CALL FAILED. TERMINATING..." + bcolors.END)
            exit(1)
