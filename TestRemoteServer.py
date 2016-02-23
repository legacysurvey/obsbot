from __future__ import print_function
import Pyro.core

#class DecamQueue(object):
class DecamQueue(Pyro.core.ObjBase):
    def execute(self, command):
        print('Executing:', command)

Pyro.core.initServer()

daemon = Pyro.core.Daemon()
uri = daemon.connect(DecamQueue(), 'CMDSRV')

# Pyro4:
#daemon = Pyro.Daemon()
#uri = daemon.register(DecamQueue)

print('Ready. Object uri =', uri)
daemon.requestLoop()
