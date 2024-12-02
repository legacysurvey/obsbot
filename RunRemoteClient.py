import os

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
