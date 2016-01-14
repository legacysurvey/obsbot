from django.db import models

class MeasuredCCD(models.Model):
    filename = models.CharField(max_length=1024)
    extension = models.CharField(max_length=32)
    expnum = models.IntegerField()
    exptime = models.FloatField()
    airmass = models.FloatField()
    racenter = models.FloatField()
    deccenter = models.FloatField()
    rabore = models.FloatField()
    decbore = models.FloatField()
    band = models.CharField(max_length=256)
    ebv = models.FloatField()
    zeropoint = models.FloatField()
    transparency = models.FloatField()
    seeing = models.FloatField()
    sky = models.FloatField()
    expfactor = models.FloatField()
