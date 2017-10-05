from __future__ import print_function
import unittest
import os
import ephem

faketime = ephem.Date('2017/10/05 01:51:52.00')
def bigben():
    return faketime


class TestMosbot2(unittest.TestCase):

    def setUp(self):
        self.testdatadir = os.path.join(os.path.dirname(__file__),
                                        'testdata')

    def test_ha(self):
        import tempfile
        from mosbot import main
        # Override the ephem.now() function to lie about the time
        oldnow = ephem.now
        ephem.now = bigben

        tempdir = tempfile.mkdtemp()
        
        fn1 = os.path.join(self.testdatadir, '2017-10-04-p3.json')

        args = ['--script', os.path.join(tempdir, 'tonight.sh')]
        args += [fn1, fn1, fn1]
        mosbot = main(cmdlineargs=args, get_mosbot=True)

        f = open(os.path.join(tempdir, 'seqnum.txt'), 'w')
        f.write('1')
        f.close()

        mosbot.Nahead = 200
        
        mosbot.update_for_image({
            'airmass': 1.019,
            'skybright': 19.372964228396995,
            'transparency': 0.94911683705518535,
            'zp': 26.348159194305978,
            'seeing': 0.9776841626480679,
            'band': 'z'})

        ephem.now = oldnow
