# -*- coding: utf-8 -*-
# Generated by Django 1.9.1 on 2016-02-16 10:01
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('obsdb', '0006_measuredccd_nmatched'),
    ]

    operations = [
        migrations.AddField(
            model_name='measuredccd',
            name='object',
            field=models.CharField(default=b'', max_length=64),
        ),
        migrations.AddField(
            model_name='measuredccd',
            name='passnumber',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='measuredccd',
            name='tileebv',
            field=models.FloatField(default=0),
        ),
        migrations.AddField(
            model_name='measuredccd',
            name='tileid',
            field=models.IntegerField(default=0),
        ),
    ]