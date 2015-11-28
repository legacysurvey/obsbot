from __future__ import print_function
from RemoteClient import RemoteClient

rc = RemoteClient()
pid = rc.execute('get_propid')
print('Got propid:', pid)

