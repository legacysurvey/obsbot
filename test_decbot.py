from __future__ import print_function
import unittest
import os
import datetime
import time

import Pyro.core
class TestQueue(Pyro.core.ObjBase):
    def __init__(self):
        super(TestQueue, self).__init__()
        self.commands = []
        self.exposures = []
    def execute(self, command):
        self.commands.append(command)
        print(type(command))
        if command.startswith('command=addexposure'):
            self.exposures.append(command)

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
        
class TestDecbot(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')
        self.server = TestServer()
        print('URL type:', type(self.server.url))

        import tempfile

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

        #self.assertEqual(len(decbot.planned_tiles), 2)
        self.assertEqual(decbot.seqnum, 2)
        self.assertEqual(len(self.server.queue.exposures), 2)

    def test_queuetime(self):
        self.server.queue.exposures = []

        from decbot import main
        args = self.jsonfiles
        args += ['--remote-server', self.server.url.address]
        args += ['--remote-port',   str(self.server.url.port)]
        args += ['--exptime', '2']
        
        print('Creating Decbot...')
        decbot = main(cmdlineargs=args, get_decbot=True)

        decbot.nom.overhead = 1.
        decbot.queueMargin = 2.
        decbot.queue_initial_exposures()
        
        self.assertEqual(decbot.seqnum, 2)
        #self.assertEqual(len(decbot.planned_tiles), 2)
        self.assertEqual(len(self.server.queue.exposures), 2)
        
        now = datetime.datetime.utcnow()
        dt = (decbot.queuetime - now).total_seconds()
        print('Time until queuetime:', dt)

        # Should be 4, minus execution time
        self.assertGreater(dt, 3)
        self.assertLess(dt, 4)

        # Should do nothing
        decbot.heartbeat()

        self.assertEqual(decbot.seqnum, 2)
        self.assertEqual(len(self.server.queue.exposures), 2)

        time.sleep(4)

        # Should queue exposure
        decbot.heartbeat()

        self.assertEqual(decbot.seqnum, 3)
        self.assertEqual(len(self.server.queue.exposures), 3)
        
        # Now the queue time should be... exp 2 + overhead 1 - margin 2 = 1
            
        now = datetime.datetime.utcnow()
        dt = (decbot.queuetime - now).total_seconds()
        print('Time until queuetime:', dt)

        # Should be 1, minus execution time
        self.assertGreater(dt, 0)
        self.assertLess(dt, 1)

if __name__ == '__main__':
    unittest.main()
