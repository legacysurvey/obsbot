from __future__ import print_function
import numpy as np
from astrometry.util.fits import fits_table
#from astrometry.util.starutil_numpy import *
import json
from camera_mosaic import ephem_observer
import ephem
from jnox import ra2hms, dec2dms
from obsbot import get_airmass

import pylab as plt

T = fits_table('obstatus/mosaic-tiles_obstatus.fits')
T.cut(T.get('pass') == 1)
T.cut(T.in_desi == 1)
T.cut(T.dec > 77.5)
#T.cut(T.dec > 80)
T.cut(T.z_done == 0)
print(len(T), 'candidate tiles')
print('RA range', T.ra.min(), T.ra.max())

#T.cut(T.ra > 155.)
#T.cut(T.dec < 81.5)

Tall = T[np.lexsort((T.ra, T.dec))]

Tall.taken = np.zeros(len(Tall), bool)

# sort by RA
I = np.argsort(T.ra)

#### HACK -- how many tiles do you need?
# take last N
I = I[-54:]
T.cut(I)
print('Cut to', len(T))
print('RA range', T.ra.min(), T.ra.max())

obs = ephem_observer()

J = json.loads(open('pass1.json').read())
print(len(J), 'tiles in pass 1 plan')

#A = json.loads(open('pass1_byhand_v2.json').read())

nexti = 0
newJ = []

for j in J:
    print()
    obs.date = ephem.Date(str(j['approx_datetime']))

    rastr  = ra2hms (j['RA' ])
    decstr = dec2dms(j['dec'])
    ephemstr = str('%s,f,%s,%s,20' % (j['object'], rastr, decstr))
    #print(ephemstr)
    etile = ephem.readdb(ephemstr)
    etile.compute(obs)

    print('Tile', j['object'], 'at', j['approx_datetime'],
          'RA,Dec', j['RA'],j['dec'])
    airmass = get_airmass(float(etile.alt))
    print('Airmass:', airmass)
    lst = obs.sidereal_time()
    print('LST', lst)
    lsthr = np.rad2deg(float(obs.sidereal_time())) / 15.
    #    print('lst', lst)


    #if lsthr < 12. and lsthr > 25./60. and lsthr < 2.+25./60.: #lsthr < 3.+37./60:  # 00:25

    #### LST range to replace
    #if lsthr < 12. and lsthr > 4.5 and lsthr < 6.25:

    time = ephem.Date(str(j['approx_datetime']))
    t0 = ephem.Date('2017-11-27 05:10:00')
    t1 = ephem.Date('2017-11-27 06:15:00')

    if time > t0 and time < t1:
        print('Over the pole')

        if False:
            # Take Arjun's script and plug in the approx_datetime
            if nexti >= len(A):
                break
            jnew = A[nexti]
            jnew.update(approx_datetime=j['approx_datetime'])
            nexti += 1
        elif False:
            # Pull from my list of tiles computed at the beginning of this
            # script
            t = T[nexti]
            nexti += 1
            jnew = dict(RA=t.ra, 
                        dec=t.dec,
                        approx_datetime=j['approx_datetime'],
                        object='MzLS_%i_z' % t.tileid,
                        filter='zd',
                        expType='object',
                        seqid='1',
                        seqnum=nexti,
                        expTime=80)
        else:
            Tall.airmass = np.zeros(len(Tall))
            Tall.ha = np.zeros(len(Tall))
            for i,t in enumerate(Tall):
                rastr  = ra2hms (t.ra)
                decstr = dec2dms(t.dec)
                obj = 'MzLS_%i_z' % t.tileid
                ephemstr = str('%s,f,%s,%s,20' % (obj, rastr, decstr))
                etile = ephem.readdb(ephemstr)
                etile.compute(obs)
                airmass = get_airmass(float(etile.alt))
                e_ra = ephem.degrees(str(t.ra))
                ha = lst - e_ra
                # over the pole?
                ha += np.pi
                Tall.airmass[i] = airmass
                Tall.ha[i] = ha
            #print('  Tile', obj, 'RA,Dec', t.ra, t.dec,
            #      'Airmass', airmass, 'HA', np.rad2deg(float(ha)), 'deg')
            I = np.flatnonzero(np.logical_not(Tall.taken) *
                               (Tall.airmass < 2.5) *
                               (Tall.ha > -80) * (Tall.ha < +80))
            print(len(I), 'tiles in airmass and HA range and not taken yet')
            for t in Tall[I]:
                print('  Tile', t.ra, t.dec, 'airmass', t.airmass, 'ha', t.ha)
            plt.clf()
            plt.plot(Tall.ra, Tall.dec, 'k.', alpha=0.5)
            K = np.flatnonzero(Tall.taken)
            plt.plot(Tall.ra[K], Tall.dec[K], 'rx', mew=2, ms=8)
            plt.scatter(Tall.ra[I], Tall.dec[I], c=Tall.airmass[I])
            plt.colorbar()
            plt.savefig('otp-%02i.png' % nexti)

            # Choose the one with largest RA
            ii = np.argmax(Tall.ra[I])
            i = I[ii]
            t = Tall[i]
            Tall.taken[i] = True
            nexti += 1
            jnew = dict(RA=t.ra, 
                        dec=t.dec,
                        approx_datetime=j['approx_datetime'],
                        object='MzLS_%i_z' % t.tileid,
                        filter='zd',
                        expType='object',
                        seqid='1',
                        seqnum=nexti,
                        expTime=80)
            
        print('tile', nexti)

        newJ.append(jnew)

        rastr  = ra2hms (jnew['RA' ])
        decstr = dec2dms(jnew['dec'])
        ephemstr = str('%s,f,%s,%s,20' % (jnew['object'], rastr, decstr))
        etile = ephem.readdb(ephemstr)
        etile.compute(obs)
        print('Tile', jnew['object'], 'RA,Dec', jnew['RA'],jnew['dec'])
        airmass = get_airmass(float(etile.alt))
        print('Airmass:', airmass)
        e_ra = ephem.degrees(str(jnew['RA']))
        #print('ephem RA:', e_ra)
        #print('in floats: RA', float(e_ra), 'LST', float(lst))
        ha = lst - e_ra
        # over the pole?
        ha += np.pi
        #if ha < -np.pi:
        #    ha += 2.*np.pi
        #print('HA', float(ha))
        print('HA:', np.rad2deg(float(ha)), 'deg')



            
open('otp.json','wb').write(json.dumps(newJ, sort_keys=True,
                                       indent=4, separators=(',', ': ')))
