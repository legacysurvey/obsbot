    
def django_setup(database_filename=None):
    import os
    from django.conf import settings
    if database_filename is None:
        database_filename = 'obsdb.sqlite3'
    exptime_database_filename = 'exptime.sqlite3'
    basedir = os.path.dirname(os.path.dirname(__file__))
    settings.configure(
        INSTALLED_APPS=['obsdb'],
        MIDDLEWARE_CLASSES=[],
        DATABASES=dict(
            default=dict(
                ENGINE='django.db.backends.sqlite3',
                NAME=os.path.join(basedir, 'obsdb', database_filename)),
            exptime=dict(
                ENGINE='django.db.backends.sqlite3',
                NAME=os.path.join(basedir, 'obsdb', exptime_database_filename)),
        ),
        ROOT_URLCONF='obsdb.urls',
        DATABASE_ROUTERS = ['obsdb.models.DatabaseRouter'],
    )

    import django
    django.setup()

    global MeasuredCCD
    from obsdb import models
    MeasuredCCD = models.MeasuredCCD

    return settings

MeasuredCCD = None
