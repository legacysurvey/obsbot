from django.db import models
import numpy as np

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

    band = models.CharField(max_length=256, default='')
    ebv = models.FloatField(default=0)
    zeropoint = models.FloatField(default=0)
    transparency = models.FloatField(default=0)
    seeing = models.FloatField(default=0)
    sky = models.FloatField(default=0)
    expfactor = models.FloatField(default=0)

    # pixel offsets vs PanSTARRS1
    dx = models.FloatField(default=0)
    dy = models.FloatField(default=0)

    # md5sum of the first hdu's data.
    # (actually, for Mosaic3 we just SUM the image pixels!)
    md5sum = models.CharField(max_length=128, default='')

    # MOSAIC bad pixel count flag set?
    bad_pixcnt = models.BooleanField(default=False)
