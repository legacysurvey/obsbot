import os
import subprocess

class RunRemoteClient(object):
    def stoprequested(self):
        cmd = 'remote-client --stop-requested'
        rtn = os.system(cmd)
        if rtn:
            print('failed to run remote-client:', rtn)

    def stopexposure(self):
        cmd = 'remote-client --stop'
        rtn = os.system(cmd)
        if rtn:
            print('failed to run remote-client:', rtn)

    def get_n_queued(self):
        cmd = 'remote-client --get-n-queued'
        cp = subprocess.run(cmd, capture_output=True, check=True, text=True)
        for line in cp.stdout.split('\n'):
            words = line.split(' ')
            if words[0] == 'n_queued':
                try:
                    n = int(words[2])
                    return n
                except:
                    pass
        print('RunRemoteClient: failed to find n_queued output.  Stdout:', cp.stdout,
              'Stderr:', cp.stderr)
        return -1

    def addexposure(self, **exp):
        pass
    def modifyexposure(self, select=None, update=None):
        pass
