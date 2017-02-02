from __future__ import print_function
from astrometry.blind.plotstuff import Plotstuff
from astrometry.util.util import Sip, anwcs, anwcs_new_sip, wcs_pv2sip_hdr
from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec
from astrometry.util.plotutils import *
import fitsio
import numpy as np

'''
avconv -r 4 -i tile-%03d.png -y mtile.mov
for x in tile*.png; do pngtopnm $x | pnmquant 256 | ppmtogif > $x.gif; done
gifsicle -d 25 -m --colors 256 tile*.gif > mtile.gif
'''

T = fits_table('obstatus/mosaic-tiles_obstatus.fits')

fn = 'k4m_160203_061211_ooi_zd_v2.fits.fz'
hdr = fitsio.read_header(fn)
obj = hdr['OBJECT']
words = obj.strip().split('_')
tileid = int(words[1])
print('Tile', tileid)
i = np.flatnonzero(T.tileid == tileid)[0]
tile = T[i]
print('Tile RA,Dec', tile.ra, tile.dec)

I,J,d = match_radec(np.array([tile.ra]), np.array([tile.dec]), T.ra, T.dec, 3.) 
Tnear = T[J]
print('Found', len(Tnear), 'tiles nearby')
# sort by distance from center
I = np.argsort(d)
Tnear.cut(I)
Tnear.rename('pass', 'passnum')

# MzLS_594295_z

ps = PlotSequence('tile', format='%03i')

wcses = []
for ext in range(1, 4+1):
    hdr = fitsio.read_header(fn, ext=ext)
    wcs = wcs_pv2sip_hdr(hdr)
    wcses.append(wcs)
    print('WCS', wcs)

PW,PH = 800,800
plot = Plotstuff(size=(PW, PH), rdw=(tile.ra, tile.dec, 2), outformat='png')
plot.color = 'verydarkblue'
plot.plot('fill')

plot.color = 'white'
plot.alpha = 0.25
plot.outline.fill = True
plot.apply_settings()

for passnum in [1, 2, 3]:

    I = np.flatnonzero(Tnear.passnum == passnum)
    for ird,(r,d) in enumerate(zip(Tnear.ra[I], Tnear.dec[I])):

        overlap = False

        allxx,allyy = [],[]

        for wcs in wcses:
            wcs.set_crval((r, d))

            W,H = wcs.get_width(), wcs.get_height()
            corners = [(1,1),(1,H),(W,H),(W,1)]
            rd = [wcs.pixelxy2radec(x,y) for x,y in corners]
            #print('rd', rd)
            ok,xx,yy = plot.wcs.radec2pixelxy(
                [ri for ri,di in rd], [di for ri,di in rd])
            allxx.append(xx)
            allyy.append(yy)
            if max(xx) < 0 or max(yy) < 0 or min(xx) > PW or min(yy) > PH:
                # No overlap; skip
                continue

            overlap = True
            #plot.outline.wcs = anwcs_new_sip(wcs)
            #plot.plot('outline')

        #if passnum == 1 and ird == 0:
        #    plot.write('pass-0.png')

        if not overlap:
            continue

        allxx,allyy = np.hstack(allxx), np.hstack(allyy)
        xlo,xhi = allxx.min(), allxx.max()
        ylo,yhi = allyy.min(), allyy.max()
        poly = [(xlo,ylo),(xlo,yhi),(xhi,yhi),(xhi,ylo)]
        print('Polygon', poly)
        plot.polygon(poly)
        plot.close_path()
        plot.fill()
        #plot.stroke()
        
        plot.write(ps.getnext())
            
    plot.write('pass-%i.png' % passnum)
            
# plot.rgb = (0.2,0.2,0.2)
# plot.plot_grid(0.1, 0.1)
# plot.color = 'gray'
# plot.plot_grid(0.5, 0.5, 0.5, 0.5)
#plot.write('plot.png')

