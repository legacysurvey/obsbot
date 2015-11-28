#
# Client for SISPI Remote CommandServer
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
            return self.cs.execute('command=%s, params = %s' % (command, parameter))

