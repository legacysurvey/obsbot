# -*- coding: utf-8 -*-
# Generated by Django 1.9.1 on 2016-02-03 23:19
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('obsdb', '0002_measuredccd_camera'),
    ]

    operations = [
        migrations.AddField(
            model_name='measuredccd',
            name='dx',
            field=models.FloatField(default=0),
        ),
        migrations.AddField(
            model_name='measuredccd',
            name='dy',
            field=models.FloatField(default=0),
        ),
        migrations.AddField(
            model_name='measuredccd',
            name='md5sum',
            field=models.CharField(default=b'', max_length=128),
        ),
    ]
