import os
import subprocess
import json

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
        cmd = ['remote-client', '--get-n-queued']
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
        cmd = ['remote-client', '--add-exposure',  json.dumps(exp)]
        rtn = os.system(cmd)
        if rtn:
            print('failed to run command', cmd, ': return val', rtn)

    def modifyexposure(self, select=None, update=None):
        cmd = ['remote-client', '--modify-exposure',  json.dumps(select), json.dumps(update)]
        rtn = os.system(cmd)
        if rtn:
            print('failed to run command', cmd, ': return val', rtn)

