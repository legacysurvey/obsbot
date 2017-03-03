from astrometry.util.fits import *
import pylab as plt
import numpy as np
import ephem
import json

ra0, ra1 = 155,175
dec0,dec1 = 15,30

def axes():
    plt.axis([ra0-1, ra1-1, dec0-1.5, dec1+1.5])
    ax = plt.axis()
    plt.xlim(ax[1],ax[0])
    plt.xticks(range(155, 180, 5))
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    

T = fits_table('obstatus/decam-tiles_obstatus.fits')

assert(np.all(T.tileid == np.arange(len(T)) + 1))

I1 = np.flatnonzero((T.ra > ra0-1) * (T.ra < ra1+1) * (T.dec > dec0-1) * (T.dec < dec1+1) * (T.get('pass') == 1))

# Tileids for passes 1,2,3 are offset by:
dt = 15872

plt.clf()
plt.plot(T.ra[I1], T.dec[I1], 'k.', alpha=0.5)
axes()
plt.savefig('/tmp/1.png')

I2 = I1 + dt
I3 = I1 + 2*dt

# plt.clf()
# plt.plot([T.ra[I1], T.ra[I2]], [T.dec[I1], T.dec[I2]], 'k-', alpha=0.5)
# plt.plot(T.ra[I2], T.dec[I2], 'b.', alpha=0.5)
# plt.plot(T.ra[I3], T.dec[I3], 'r.', alpha=0.5)
# plt.savefig('/tmp/2.png')

P1 = T[I1]
P2 = T[I2]
P3 = T[I3]

plt.clf()
plt.plot(T.ra, T.dec, 'k.', alpha=0.5)
plt.plot(T.ra[T.r_done == 1], T.dec[T.r_done == 1] + 0.2, 'r.', alpha=0.5)
plt.plot(T.ra[T.g_done == 1], T.dec[T.g_done == 1] - 0.2, 'g.', alpha=0.5)
axes()
plt.title('DONE tiles')
plt.savefig('/tmp/3.png')

Ir = np.flatnonzero((P1.r_done == 0) * (P2.r_done == 0) * (P3.r_done == 0) *
                    (P1.ra > ra0) * (P1.ra < ra1) *
                    (P1.dec > dec0) * (P1.dec < dec1))
Ig = np.flatnonzero((P1.g_done == 0) * (P2.g_done == 0) * (P3.g_done == 0) *
                    (P1.ra > ra0) * (P1.ra < ra1) *
                    (P1.dec > dec0) * (P1.dec < dec1))

print('Need', len(Ir), 'r and', len(Ig), 'g')

plt.clf()
plt.plot(T.ra[I1], T.dec[I1], 'k.', alpha=0.5)
plt.plot(P1.ra[Ir], P1.dec[Ir] + 0.2, 'r.', alpha=0.5)
plt.plot(P1.ra[Ig], P1.dec[Ig] - 0.2, 'g.', alpha=0.5)
axes()
plt.title('Needed tiles')
plt.savefig('/tmp/4.png')



PP = P1[np.hstack([Ir, Ig])]
PP.band = np.array(['r'] * len(Ir) + ['g'] * len(Ig))

#K = np.argsort(PP.ra)

K1 = np.flatnonzero(PP.dec < 18)
K2 = np.flatnonzero((PP.dec >= 18) * (PP.dec <= 25.))
K3 = np.flatnonzero(PP.dec > 25)
K1 = K1[np.argsort( PP.ra[K1])]
K2 = K2[np.argsort(-PP.ra[K2])]
K3 = K3[np.argsort( PP.ra[K3])]
K = np.hstack([K1,K2,K3])

PP.cut(K)

exptimes = dict(g=56., r=40.)

date = ephem.Date('2017/03/02 03:45:00')

overhead = 30.

plan = []
for i,tile in enumerate(PP):
    plan.append(dict(seqid='1', seqnum=i+1, seqtot=len(PP),
                     expType='object',
                     object='DECaLS_%i_%s' % (tile.tileid, tile.band),
                     expTime=exptimes[tile.band],
                     filter=tile.band,
                     approx_datetime=str(ephem.Date(date)),
                     RA=tile.ra,
                     dec=tile.dec))
    date += (exptimes[tile.band] + overhead) / (24.*3600.)

f = open('eboss.json', 'w')
f.write(json.dumps(plan).replace(',',',\n').replace('{','\n{\n').replace('}','\n}\n'))
f.close()

