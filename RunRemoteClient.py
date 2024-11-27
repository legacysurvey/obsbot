import os

class RunRemoteClient(object):
    def stopexposure(self):
        cmd = 'remote-client --stop'
        rtn = os.system(cmd)
        if rtn:
            print('failed to run remote-client:', rtn)


