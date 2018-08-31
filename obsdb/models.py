from django.db import models
import numpy as np

class ComputedExptime(models.Model):

    def __str__(self):
        return ('ComputedExptime: start %.4f, seq %i' % (self.starttime, self.seqnum) +
                ', tile %i, pass %i, band %s' % (self.tileid, self.passnumber, self.band) +
                ', adj %.3f, exptime %.1f' % (self.adjfactor, self.exptime))

    starttime = models.FloatField()
    seqnum = models.IntegerField()

    tileid = models.IntegerField(default=0)
    passnumber = models.IntegerField(default=0)
    band = models.CharField(max_length=16, default='')
    airmass = models.FloatField(default=0)
    ebv = models.FloatField(default=0)

    # properties of the just-measured image being used to determine exptime
    meas_band = models.CharField(max_length=16, default='')
    zeropoint = models.FloatField(default=0)
    transparency = models.FloatField(default=0)
    seeing = models.FloatField(default=0)
    sky = models.FloatField(default=0)
    expfactor = models.FloatField(default=0)

    # other passes -- via ForeignKey many-to-one

    adjfactor = models.FloatField(default=0)

    exptime_unclipped = models.FloatField(default=0)
    exptime_clipped = models.FloatField(default=0)
    exptime_satclipped = models.FloatField(default=0)
    exptime = models.FloatField(default=0)

class OtherPasses(models.Model):
    exposure = models.ForeignKey(ComputedExptime, on_delete=models.DO_NOTHING)

    tileid = models.IntegerField(default=0)
    passnumber = models.IntegerField(default=0)
    depth = models.FloatField(default=0)


class DatabaseRouter(object):
    def db_for_read(self, model, **hints):
        # print('DatabaseRouter.db_for_read:', model)
        if model in [ComputedExptime, OtherPasses]:
            # print('Returned "exptime"')
            return 'exptime'
        return None
    db_for_write = db_for_read


class MeasuredCCD(models.Model):
    camera = models.CharField(
        max_length=32,
        choices=(('decam','DECam'),
                 ('mosaic3','Mosaic3')),
        default='decam')
    filename = models.CharField(max_length=1024)
    extension = models.CharField(max_length=32)
    expnum = models.IntegerField(default=0)
    exptime = models.FloatField(default=0)
    mjd_obs = models.FloatField(default=0)
    airmass = models.FloatField(default=0)
    # Approximate (header) RA,Dec center of CCD
    racenter = models.FloatField(default=0)
    deccenter = models.FloatField(default=0)
    # Approximate (header) RA,Dec center of camera
    rabore = models.FloatField(default=0)
    decbore = models.FloatField(default=0)

    obstype = models.CharField(max_length=64, default='object')
    
    object = models.CharField(max_length=64, default='')
    tileid = models.IntegerField(default=0)
    passnumber = models.IntegerField(default=0)
    tileebv = models.FloatField(default=0)

    band = models.CharField(max_length=256, default='')
    ebv = models.FloatField(default=0)
    zeropoint = models.FloatField(default=0)
    transparency = models.FloatField(default=0)
    seeing = models.FloatField(default=0)
    sky = models.FloatField(default=0)
    expfactor = models.FloatField(default=0)

    # pixel offsets vs Pan-STARRS1
    dx = models.FloatField(default=0)
    dy = models.FloatField(default=0)

    # WCS
    # cd1_1 = models.FloatField(default=0)
    # cd1_2 = models.FloatField(default=0)
    # cd2_1 = models.FloatField(default=0)
    # cd2_2 = models.FloatField(default=0)
    
    # number of stars matched to Pan-STARRS1
    nmatched = models.IntegerField(default=-1)

    # md5sum of the first hdu's data.
    # (actually, for Mosaic3 we just SUM the image pixels!)
    md5sum = models.CharField(max_length=128, default='')

    # MOSAIC bad pixel count flag set?
    bad_pixcnt = models.BooleanField(default=False)

    # MOSAIC -- read time
    readtime = models.FloatField(default=0)

    affine_dx = models.FloatField(default=0)
    affine_dxx = models.FloatField(default=0)
    affine_dxy = models.FloatField(default=0)
    affine_dy = models.FloatField(default=0)
    affine_dyx = models.FloatField(default=0)
    affine_dyy = models.FloatField(default=0)
    affine_x0 = models.FloatField(default=0)
    affine_y0 = models.FloatField(default=0)




