import subprocess
import os

def call(cmd):
    # output = subprocess.check_output(['ls', '-1'])
    # print 'Have %d bytes in output' % len(output)
    print(">>> {}".format(cmd))
    os.system(cmd)