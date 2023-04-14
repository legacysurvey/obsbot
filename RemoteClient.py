from __future__ import print_function
#
# Client for SISPI Remote CommandServer
# Code from Klaus Honscheid, 2015-11-06, [decam-chatter 1546]
#
import Pyro.core

class RemoteClient():
    def __init__(self, cs_host='system1', cs_port=15555, cs_name = 'CMDSRV'):
        """ Initialize remote client class """
        self.host = str(cs_host)
        self.port = int(cs_port)
        self.name = str(cs_name)
        self.uri = "PYROLOC://%s:%d/%s" % (self.host, self.port, self.name)

        # Connect to command server. Exceptions to be handled by caller
        self.cs = Pyro.core.getProxyForURI(self.uri)
       
    def execute(self, command, parameter=None):
        """ execute remote command """
        if parameter == None:
            return self.cs.execute('command=%s' % command)
        else:
            cmd = 'command=%s, params = %s' % (command, parameter)
            print('Sending execute("%s")' % cmd)

            return self.cs.execute('command=%s, params = %s' % (command, parameter))

    # Added by Dustin -- expose specific functions Klaus has provided for us
    def get_propid(self):
        return self.execute('get_propid')

    def stopexposure(self):
        return self.execute('stopexposure')

    def get_n_queued(self):
        n = self.execute('get_nqueue')
        return int(n)

    def clear_queue(self):
        return self.execute('clear_queue')

    def addexposure(self, exptime=10., exptype='object', filter='r',
                    object=None, ra=0., dec=0., verbose=False, propid=None):
        import json
        kw = dict(expType=exptype, object=object)
        if propid is not None:
            kw.update(propid=propid)
        if exptype == 'object':
            if object is None:
                object = 'Object'
            kw.update(expTime=exptime, filter=filter, RA=ra, dec=dec)
            paramstr = json.dumps(kw)
                                       
        elif exptype == 'dark':
            if object is None:
                object = 'dark'
            kw.update(expTime=exptime)
            paramstr = json.dumps(kw)

        elif exptype == 'dome flat':
            if object is None:
                object = 'flat'
            kw.update(expTime=exptime, filter=filter)
            paramstr = json.dumps(kw)

        elif exptype == 'zero':
            if object is None:
                object = 'zero'
            paramstr = json.dumps(kw)
            #paramstr = 'exptype=%s' % exptype

        if verbose:
            print("Calling addexposure, parameter='%s'" % paramstr)

        return self.execute('addexposure', parameter=paramstr)



if __name__ == '__main__':
    # rc = RemoteClient(cs_host='10.10.168.162', cs_port=7767)
    # rc.addexposure()
    # rc.addexposure(exptype='dark')
    # rc.addexposure(exptype='zero')
    # rc.addexposure(exptype='dome flat', filter='z', exptime=50)

    rc = RemoteClient()
    #res = rc.addexposure(exptype='zero')
    #res = rc.addexposure(ra=42, dec=12, exptime=17., filter='g')

    res = rc.get_propid()
    print('Got propid:', res)

    res = rc.get_n_queued()
    print('Got N queued:', res)

    res = rc.addexposure(ra=80., dec=-10., verbose=True,
                         propid='2014B-0404',
                         #propid='2023A-140687',
                         object='DECaLS_5774_z')
    print('Addexposure: got', res)
