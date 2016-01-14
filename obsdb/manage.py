#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    from django.conf import settings
    # basedir = os.path.dirname(os.path.dirname(__file__))
    # settings.configure(INSTALLED_APPS=['obsdb'],
    #                    MIDDLEWARE_CLASSES=[],
    #                    DATABASES=dict(default=dict(
    #                        ENGINE='django.db.backends.sqlite3',
    #                        NAME=os.path.join(basedir,'obsdb','obsdb.sqlite3')))
    #     )

    import obsdb
    obsdb.django_setup()
    
    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
