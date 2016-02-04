# -*- coding: utf-8 -*-
# Generated by Django 1.9.1 on 2016-01-20 07:40
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('obsdb', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='measuredccd',
            name='camera',
            field=models.CharField(choices=[(b'decam', b'DECam'), (b'mosaic3', b'Mosaic3')], default=b'decam', max_length=32),
        ),
    ]