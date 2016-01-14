from django.db import models
import numpy as np

class MeasuredCCD(models.Model):
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
