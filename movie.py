from __future__ import print_function
import os
import datetime

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import json
import ephem

from astrometry.util.starutil_numpy import degrees_between

def main():
    import optparse
    import sys
    
    parser = optparse.OptionParser(usage='%prog <json>')
    parser.add_option('--base', default='plan', help='Plot base filename')
    parser.add_option('-t', '--obstatus', help='Show already-observed tiles?')
    parser.add_option('--bands', help='Plot only already-observed tiles in the given bands', default='g,r,z')
    parser.add_option('--sgc', action='store_true', help='Center on SGC?')

    parser.add_option('--ralo',  type=float, default=None)
    parser.add_option('--rahi',  type=float, default=None)
    parser.add_option('--declo', type=float, default=None)
    parser.add_option('--dechi', type=float, default=None)
    parser.add_option('--scaled', action='store_true', default=False,
                      help='Scale plot so that 1 deg RA = 1 deg Dec (no COS term)')
    parser.add_option('--wide', action='store_true', default=False,
                      help='Make wider plots?')

    parser.add_option('--also', action='append', default=[],
                      help='Also plot the plan from the given filename.')
    
    parser.add_option('--mosaic', action='store_true', help='Set defaults for Mosaic survey')
    
    parser.add_option('--start-time', help='Start time for this plan, HH:MM:SS UTC.  Default: 12-degree twilight tonight.')
    parser.add_option('--start-date', help='Start date for this plan, YYYY-MM-DD UTC.')
    parser.add_option('--second-half', action='store_true', help='This plan starts at the start of the second half-night.')
    
    parser.add_option('--skip', type=int, default=1, help='Write every Nth plot only')
    
    opt,args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit(-1)

    if opt.mosaic:
        dd = dict(ralo=0, rahi=360, declo=20, dechi=70)
    else:
        dd = dict(ralo=0, rahi=360, declo=-10, dechi=35)
    for k in dd.keys():
        if getattr(opt, k, None) is None:
            setattr(opt, k, dd[k])
        
    start_date_specified = (opt.start_date is not None)
        
    if opt.start_date is None:
        # Get date at start of night, where we define a new day as
        # starting at noon UTC.
        now = datetime.datetime.utcnow()
        # noon
        nightstart = now - datetime.timedelta(0, 12 * 3600)
        d = nightstart.date()
        opt.start_date = '%04i-%02i-%02i' % (d.year, d.month, d.day)
        print('Set start date to', opt.start_date)
        
    obs = ephem.Observer()
    if opt.mosaic:
        obs.lon = '111.6'
        obs.lat = '31.9633'
        obs.elev = 2064.
    else:
        obs.lon = '-70.806525'
        obs.lat = '-30.169661'
        obs.elev = 2207.0 # meters

    obs.temp = 10.0 # deg celsius; average temp for August
    obs.pressure = 780.0 # mbar

    ### HACK
    obs.date = ephem.Date(opt.start_date + ' 8:00:00')
    #print('Obs date:', obs.date)
    daystart = obs.date
    obs.horizon = -ephem.degrees('12:00:00.0')
    sun = ephem.Sun()
    eve_twi  = obs.next_setting(sun)
    obs.date = eve_twi
    morn_twi = obs.next_rising(sun)
    print('Evening twilight:', eve_twi)
    print('Morning twilight:', morn_twi)
    assert(morn_twi > eve_twi)
    obs.horizon = 0.
    print('Eve twi:', eve_twi, 'Morning:', morn_twi)
    
    if opt.second_half:
        # Set start-time to the midpoint between 12-degree twilights.
        obs.date = ephem.Date((eve_twi + morn_twi) / 2.)
        print('Second half starts at', obs.date)
        
    elif opt.start_time is None:
        # 12-degree twilight on start_date
        obs.date = eve_twi

    else:
        obs.date = ephem.Date(opt.start_date + ' ' + opt.start_time)
        if not start_date_specified and obs.date < daystart:
            # If --start-date is, eg, 2am, assume it's during the night starting on daystart.
            obs.date += 1.
        print('Start date:', obs.date)

    jfn = args[0]
    print('Reading JSON file', jfn)
    J = json.loads(open(jfn,'rb').read())
    print(len(J), 'entries')

    Jalso = [json.loads(open(fn,'rb').read()) for fn in opt.also]

    tiles = None
    if opt.obstatus is not None:
        from astrometry.util.fits import fits_table
        
        tiles = fits_table(opt.obstatus)
        print('Read', len(tiles), 'tiles')
        tiles = tiles[(tiles.in_des == 0) * np.logical_or(
            (tiles.in_sdss == 1),
            (tiles.in_sdss == 0) * (tiles.in_desi == 1))]
        print(len(tiles), 'in footprint')
    
    fcmap = dict(g='g',r='r',z='m', zd='m')
    ddecmap = dict(g=-0.2, r=0, z=0.2, zd=0.2)
    
    ras = np.array([j['RA'] for j in J])
    decs = np.array([j['dec'] for j in J])
    filts = np.array([j['filter'] for j in J])
    exptime = np.array([j['expTime'] for j in J])
    fieldname = [j['object'] for j in J]
    passnum = np.zeros(len(J), int)

    #
    # Observe the horribleness of --sgc:
    # We LIE to matplotlib, being careful to use transform_ra() to go from
    # real RA to plotted x coordinate; we add 180 and mod by 360.  Ick!
    # We then explicitly set the tick marks to hide the body... ahem, cover
    # our tracks.
    #
    SGC_DRA = 180.
    if opt.sgc:
        transform_ra = lambda x: (x + SGC_DRA) % 360.
    else:
        transform_ra = lambda x: x
    
    filtcc = np.array([fcmap[f] for f in filts])
    ddecs = np.array([ddecmap[f] for f in filts])

    # passmap = { 1: dict(marker='.'),
    #             2: dict(marker='o', mfc='none'),
    #             3: dict(marker='x') }
    passmap = { 1: dict(marker='.'),
                2: dict(marker='.'),
                3: dict(marker='.'), }

    opt.bands = opt.bands.split(',')
    if len(opt.bands) == 1:
        filtddec = {'g':0, 'r':0, 'z':0}
    else:
        ddec = 0.4
        filtddec = { 'g': -ddec, 'r': 0, 'z': ddec }
    
    seqmap = ['r','y','g','b','m']
    #seqcc = np.array([seqmap[s % len(seqmap)] for s in seqnum])
    #seqcc = np.array([seqmap[s % len(seqmap)] for s in seqid])
    
    if opt.wide:
        plt.figure(figsize=(12,8))
    plt.subplots_adjust(left=0.1, right=0.95)

    plt.clf()
    plt.plot(transform_ra(ras), decs, 'r.')
    plt.axis('scaled')
    ax = [opt.rahi, opt.ralo, opt.declo, opt.dechi]
    plt.axis(ax)

    moon = ephem.Moon()

    # Get times when exposures should occur.
    times = []
    LSTs = []

    # If the JSON files include estimated times, use those
    if 'approx_datetime' in J[0]:
        for j in J:
            obs.date = ephem.Date(str(j['approx_datetime']))
            times.append(ephem.Date(obs.date))
            LSTs.append(np.rad2deg(float(obs.sidereal_time())))
    else:
        # Predict overheads
        lastra,lastdec = None,None
        for i in range(len(J)):
            print('Exposure', i, 'should start at', str(obs.date))
            times.append(ephem.Date(obs.date))
            LSTs.append(np.rad2deg(float(obs.sidereal_time())))
            overhead = 30.
            if lastra is not None:
                slew = degrees_between(lastra, lastdec, ras[i], decs[i])
                lastra  = ras [i]
                lastdec = decs[i]
                # Add 3 seconds per degree for slews longer than 2 degrees
                overhead += np.maximum(0, slew - 2.) * 3.
            # Add overhead
            print('Adding', exptime[i], 'seconds exptime plus',
                  overhead, 'seconds overhead')
            obs.date += (exptime[i] + overhead) / (24 * 3600.)

    also = []
    for Ja in Jalso:
        atimes = np.array([ephem.Date(str(j['approx_datetime'])) for j in Ja])
        aras = np.array([j['RA'] for j in Ja])
        adecs = np.array([j['dec'] for j in Ja])
        afilts = np.array([j['filter'] for j in Ja])
        aexptime = np.array([j['expTime'] for j in Ja])
        afieldname = [j['object'] for j in Ja]
        also.append((atimes, aras, adecs, afilts, aexptime, afieldname))
            
    # Try to get the pass number via parsing the field name to get tile id
    # and looking up the pass number in the tiles table.
    if tiles is not None:
        for i,f in enumerate(fieldname):
            # "DECaLS_7884_z"
            parts = f.split('_')
            if len(parts) < 3:
                continue
            if parts[0] != 'DECaLS':
                continue
            tileid = int(parts[1])
            I = np.flatnonzero(tileid == tiles.tileid)
            if len(I) != 1:
                continue
            pa = tiles.get('pass')[I[0]]
            print('Field', f, 'tileid', tileid, 'pass', pa)
            passnum[i] = pa

    for i in reversed(range(0,len(J),opt.skip)):

        print('Exposure', i, 'of', len(J))
        plt.clf()

        if tiles is not None:
            todo = np.ones(len(tiles), bool)
            for filt in opt.bands:
                for p in [1,2,3]:
                    I = ((tiles.get('pass') == p) *
                         (tiles.get('%s_done' % filt) == 1))
                    todo[I] = False
                    #print sum(I), 'tiles done in', filt, 'pass', p
                    plt.plot(transform_ra(tiles.ra[I]), tiles.dec[I] + filtddec[filt],
                             linestyle='none',
                             color=fcmap[filt], alpha=0.25, mec=fcmap[filt],
                             zorder=10,
                             **passmap[p])
            #print sum(todo), 'tiles to-do'
            aa = 0.1
            if opt.mosaic:
                aa = 0.03
            plt.plot(transform_ra(tiles.ra[todo]), tiles.dec[todo], '.',
                     color='k', alpha=aa, zorder=15)

        rr = ras[:i+1]
        dd = decs[:i+1] + ddecs[:i+1]
        plt.scatter(transform_ra(rr), dd, c=filtcc[:i+1], s=40, zorder=50,
                    edgecolors='k')
        plt.plot(transform_ra(rr), dd, 'k-', alpha=0.5, zorder=40)

        plt.text(transform_ra(rr[i]-5), dd[i], fieldname[i],
                 bbox=dict(facecolor='w', alpha=0.8, edgecolor='none'),
                 zorder=60)

        for (atimes, aras, adecs, afilts, aexptime, afieldname) in also:
            I = np.flatnonzero(atimes <= times[i])
            afiltcc = np.array([fcmap[f]  for f in afilts])
            addecs = np.array([ddecmap[f] for f in afilts])
            rr = aras[I]
            dd = adecs[I] + ddecs[I]
            plt.scatter(transform_ra(rr), dd, c=filtcc[:i+1], s=40, zorder=50,
                        edgecolors='k', alpha=0.5)
            plt.plot(transform_ra(rr), dd, 'k-', alpha=0.25, zorder=40)
            ii = I[-1]
            plt.text(transform_ra(rr[-1]-5), dd[-1], afieldname[ii],
                     bbox=dict(facecolor='w', alpha=0.8, edgecolor='none'),
                     zorder=50)

        
        print('time:', times[i])
        obs.date = times[i]
        moon.compute(obs)
        moonra  = np.rad2deg(moon.ra)
        moondec = np.rad2deg(moon.dec)
        print('Moon RA,Dec', moonra, moondec)
        plt.plot(transform_ra(moonra), moondec, 'o', ms=20, mec=(1,0.6,0), mew=5, mfc='none', zorder=40)
        #plt.plot(moonra, moondec, 'o', ms=20, mec='k', mew=1)

        # Plot airmass contours
        dd = np.linspace(ax[2], ax[3], 20)
        rr = np.linspace(ax[0], ax[1], 21)
        airmass = np.zeros((len(dd), len(rr)))
        moonsep = np.zeros((len(dd), len(rr)))
        for ii,ra in enumerate(rr):
            for jj,dec in enumerate(dd):
                ra_str  = ConvertRA (np.array([ra ]))[0]
                dec_str = ConvertDec(np.array([dec]))[0]
                key = '%i,f,%s,%s,20' % (0, ra_str, dec_str)
                grid = ephem.readdb(key)
                grid.compute(obs)
                am = GetAirmass(float(grid.alt))
                airmass[jj,ii] = am

                ms = np.rad2deg(ephem.separation((moon.az, moon.alt),
                                                 (grid.az, grid.alt)))
                moonsep[jj,ii] = ms
                    
        levels = np.append(np.arange(1.0, 2.5, 0.1), [2.5, 3.0, 4.0])
        darkblue = (0.03, 0.19, 0.42, 0.5)
        if opt.sgc:
            # Plot the contours in two parts... if SGC_DRA != an rr grid point, this
            # may be ugly...
            I1 = np.flatnonzero(rr <= SGC_DRA)
            trr1 = transform_ra(rr[I1])
            trr1 += (trr1 < 180) * 360.
            I2 = np.flatnonzero(rr >= SGC_DRA)
            trr2 = transform_ra(rr[I2])

            for trr,am in [(trr1, airmass[:,I1]), (trr2, airmass[:,I2])]:
                plt.contourf(trr, dd, am, levels,
                             cmap='Blues', alpha=0.2, vmin=1.0, vmax=2.2)
                con = plt.contour(trr, dd, am, levels, colors=[darkblue], alpha=0.1)
                plt.clabel(con, inline=1, fontsize=10, fmt='%.1f',
                           use_clabeltext=True, alpha=0.25)
                plt.contour(trr, dd, am, [2.0, 2.5], colors=[darkblue], linewidths=[2])
            
        else:
            plt.contourf(transform_ra(rr), dd, airmass, levels, cmap='Blues', alpha=0.2,
                         vmin=1.0, vmax=2.2)
            con = plt.contour(transform_ra(rr), dd, airmass, levels, colors=[darkblue], alpha=0.1)
            plt.clabel(con, inline=1, fontsize=10, fmt='%.1f',
                   use_clabeltext=True, alpha=0.25)
            plt.contour(transform_ra(rr), dd, airmass, [2.0, 2.5], colors=[darkblue], linewidths=[2])

        #levels = np.array([10,20,30,40,50,60,70,80])
        levels = np.array([40,50,60])
        #plt.contourf(rr, dd, moonsep, levels, cmap='Blues', alpha=0.2,
        #            vmin=1.0, vmax=2.2)
        #darkblue = (0.03, 0.19, 0.42, 0.5)
        con = plt.contour(transform_ra(rr), dd, moonsep, levels, colors='r')
        plt.clabel(con, inline=1, fontsize=10, fmt='%i',
                   use_clabeltext=True, alpha=0.25)

        LST = LSTs[i]
        plt.axvline(transform_ra(LST), color='0.5')
            
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        tt = ('%s: (%.1f,%.1f), pass %i, UT: %s, %i sec' %
              (fieldname[i], ras[i], decs[i],passnum[i], times[i], exptime[i]))
        for (atimes, aras, adecs, afilts, aexptime, afieldname) in also:
            I = np.flatnonzero(atimes <= times[i])
            ii = I[-1]
            tt += ('\n%s, (%.1f, %.1f), %i sec' %
                   (afieldname[ii], aras[ii], adecs[ii], aexptime[ii]))

        plt.title(tt)

        fn = '%s-%03i.png' % (opt.base, i)

        tt = np.arange(0, 361, 60)
        if opt.rahi - opt.ralo <= 120:
            tt = np.arange(0, 361, 20)
        if opt.sgc:
            plt.xticks(transform_ra(tt), ['%i' % t for t in tt])
        else:
            plt.xticks(tt)
               
        plt.axis(ax)
        if opt.scaled:
            plt.axis('scaled')
            plt.axis(ax)
        plt.savefig(os.path.join(os.path.dirname(args[0]),fn))
        print('Wrote', fn)
        
    print()
    cmd = 'avconv -r 4 -i %s-%%03d.png -y %s.mov' % (opt.base, opt.base)
    print(cmd)
    os.system(cmd)

        
#### From nightlystrategy.py:

def GetAirmass(al):
    if (al < 0.07):
        al = 0.07
    secz = 1.0/np.sin(al)
    seczm1 = secz-1.0
    airm = secz-0.0018167*seczm1-0.002875*seczm1**2-0.0008083*seczm1**3
    return airm

def ConvertRA(raval):
    hours = np.zeros_like(raval)
    minutes = np.zeros_like(raval)
    seconds = np.zeros_like(raval)
    
    hours = (raval/360.0)*24.0
    minutes = (hours-np.floor(hours))*60.0
    seconds = (minutes-np.floor(minutes))*60.0
    
    
    stringra = []
    for k in range(0,raval.size):
        #print hours[k],minutes[k], seconds[k]
        stringra.append("%02d:%02d:%04.1f" % (hours[k], minutes[k], seconds[k]))
    
    stringra = np.array(stringra)
    return stringra


#######################################################
def ConvertDec(decval):
    sdd = np.zeros_like(decval)
    minutes = np.zeros_like(decval)
    seconds = np.zeros_like(decval)
    
    sdd = decval
    pos_sdd = np.fabs(sdd)
    minutes = (pos_sdd-np.floor(pos_sdd))*60.0
    seconds = (minutes-np.floor(minutes))*60.0
    
    stringdec = []
    for k in range(0,decval.size):
        #print sdd[k],minutes[k], seconds[k]
        stringdec.append("%02d:%02d:%02d" % (sdd[k], minutes[k], seconds[k]))
    
    stringdec = np.array(stringdec)
    return stringdec


if __name__ == '__main__':
    main()
