from __future__ import print_function
import Pyro.core

from test_decbot import TestQueue
import time

def run_queue(q):
    while True:
        e = q.pop()
        if e is None:
            print('Q empty')
            time.sleep(1)
            continue
        exptime = e['expTime']
        print('Q: Exposing:', e)
        time.sleep(exptime)
        print('Q: Overhead')
        time.sleep(q.overhead)
        print('Q: Done:', e)

Pyro.core.initServer()

testq = TestQueue()
testq.overhead = 5.

daemon = Pyro.core.Daemon()
uri = daemon.connect(testq, 'CMDSRV')

import threading
qthread = threading.Thread(target=run_queue, args=(testq,))
qthread.daemon = True
qthread.start()

print('Ready. Object uri =', uri)
daemon.requestLoop()
