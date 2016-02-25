#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    basedir = os.path.dirname(os.path.dirname(__file__))
    sys.path.append(basedir)
    
    import obsdb
    print 'obsdb:', dir(obsdb)
    obsdb.django_setup()

    from django.core.management import execute_from_command_line
    execute_from_command_line(sys.argv)
