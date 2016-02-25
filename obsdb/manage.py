#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    basedir = os.path.dirname(os.path.dirname(__file__))
    #print 'basedir:', basedir
    sys.path.append(basedir)

    from django.conf import settings
    
    import obsdb
    print 'obsdb:', dir(obsdb)
    obsdb.django_setup()

    from django.core.management import execute_from_command_line
    execute_from_command_line(sys.argv)
