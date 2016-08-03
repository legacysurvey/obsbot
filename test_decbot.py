from __future__ import print_function
import unittest
import os
import datetime
import time

# module is named 'queue' in python3
from Queue import Queue, Empty
import json

import Pyro.core

class TestQueue(Pyro.core.ObjBase):

    def __init__(self):
        super(TestQueue, self).__init__()
        self.overhead = 30.
        self.reset()

    def execute(self, command):
        self.commands.append(command)
        if command.startswith('command=addexposure'):
            # command=addexposure, params = {"filter": "z", "dec": 19.202,
            # "RA": 114.479, "expType": "object", "object": "DECaLS_23342_z",
            # "expTime": 118}
            params = command.replace('command=addexposure, params = ', '')
            params = json.loads(params)
            print('Exposure:', params)
            self.exposures.put(params)
            print('Queue now has', self.exposures.qsize(), 'entries')
            return 'SUCCESS'
        elif command.startswith('command=get_nqueue'):
            sz = self.nqueued()
            print('get_nqueue() ->', sz)
            return sz
        else:
            print('Command: "%s"' % command)

    def nqueued(self):
        return self.exposures.qsize()

    def pop(self):
        try:
            return self.exposures.get_nowait()
        except Empty:
            return None

    def reset(self):
        self.commands = []
        self.exposures = Queue()

import threading
class TestServer(object):
    def __init__(self):
        Pyro.core.initServer()
        self.queue = TestQueue()
        self.daemon = Pyro.core.Daemon()
        self.url = self.daemon.connect(self.queue, 'CMDSRV')
        print('Daemon URL:', self.url)
        self.thread = threading.Thread(target=self)
        self.thread.daemon = True
        self.thread.start()

    def __call__(self):
        print('TestServer Running at URL', self.url)
        self.daemon.requestLoop()

    def nqueued(self):
        return self.queue.nqueued()

class TestDecbot(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')
        self.server = TestServer()
        print('URL type:', type(self.server.url))

        fn1 = os.path.join(self.testdatadir, 'decals-pass1.json')
        fn2 = os.path.join(self.testdatadir, 'decals-pass2.json')
        fn3 = os.path.join(self.testdatadir, 'decals-pass3.json')

        tmpfn1 = fn1 + '.tmp'
        tmpfn2 = fn2 + '.tmp'
        tmpfn3 = fn3 + '.tmp'

        for fn,tmpfn in ((fn1,tmpfn1),(fn2,tmpfn2),(fn3,tmpfn3)):
            import json
            import ephem

            J = json.loads(open(fn, 'r').read())
            t0 = ephem.Date(str(J[0]['approx_datetime']))
            print('First exposure:', t0)
            now = ephem.now()
            print('Now:', now)
            for j in J:
                tnew = now + ephem.Date(str(j['approx_datetime'])) - t0
                tnew = str(ephem.Date(tnew))
                j['approx_datetime'] = tnew
                #print('Updated datetime to', tnew)

            f = open(tmpfn, 'w')
            json.dump(J, f, sort_keys=True, indent=4, separators=(',', ': '))
            f.close()

        self.jsonfiles = [tmpfn1, tmpfn2, tmpfn3]

    def test_new_file(self):
        from decbot import main
        args = self.jsonfiles
        args += ['--remote-server', self.server.url.address]
        args += ['--remote-port',   str(self.server.url.port)]

        decbot = main(cmdlineargs=args, get_decbot=True)
        decbot.queue_initial_exposures()

        fn = os.path.join(self.testdatadir, 'decam-00488199-n4.fits.fz')
        decbot.process_file(fn)

        self.assertEqual(len(decbot.queued_tiles), 2)
        self.assertEqual(self.server.nqueued(), 2)

    def test_nqueue(self):
        self.server.queue.reset()

        from decbot import main
        args = self.jsonfiles
        args += ['--remote-server', self.server.url.address]
        args += ['--remote-port',   str(self.server.url.port)]

        print('Creating Decbot...')
        decbot = main(cmdlineargs=args, get_decbot=True)

        decbot.queue_initial_exposures()

        self.assertEqual(len(decbot.queued_tiles), 2)
        self.assertEqual(self.server.nqueued(), 2)

        # Should do nothing
        decbot.heartbeat()

        self.assertEqual(len(decbot.queued_tiles), 2)
        self.assertEqual(self.server.nqueued(), 2)

        self.server.queue.pop()
        self.assertEqual(self.server.nqueued(), 1)

        # NOW it should queue exposure.
        decbot.heartbeat()

        self.assertEqual(len(decbot.queued_tiles), 3)
        self.assertEqual(self.server.nqueued(), 2)

    def test_ignore_missing_dir(self):
        self.server.queue.reset()

        from decbot import main
        args = self.jsonfiles
        args += ['--remote-server', self.server.url.address]
        args += ['--remote-port',   str(self.server.url.port)]

        import tempfile
        dirnm = tempfile.mkdtemp()
        os.rmdir(dirnm)
        args += ['--rawdata', dirnm]

        print('Creating Decbot...')
        decbot = main(cmdlineargs=args, get_decbot=True)

        decbot.queue_initial_exposures()

        self.assertEqual(len(decbot.queued_tiles), 2)
        self.assertEqual(self.server.nqueued(), 2)

        decbot.heartbeat()

        os.mkdir(dirnm)

        decbot.heartbeat()

        fn = os.path.join(self.testdatadir, 'decam-00488199-n4.fits.fz')
        destfn = os.path.join(dirnm, 'new.fits.fz')
        os.symlink(fn, destfn)

        self.server.queue.pop()
        self.assertEqual(self.server.nqueued(), 1)

        decbot.heartbeat()

        self.assertEqual(len(decbot.queued_tiles), 3)
        self.assertEqual(self.server.nqueued(), 2)

if __name__ == '__main__':
    unittest.main()
