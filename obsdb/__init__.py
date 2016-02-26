
def django_setup(database_filename=None):
    import os
    from django.conf import settings
    if database_filename is None:
        database_filename = 'obsdb.sqlite3'
    basedir = os.path.dirname(os.path.dirname(__file__))
    settings.configure(INSTALLED_APPS=['obsdb'],
                       MIDDLEWARE_CLASSES=[],
                       DATABASES=dict(default=dict(
                           ENGINE='django.db.backends.sqlite3',
                           NAME=os.path.join(basedir, 'obsdb',
                                             database_filename))
                      ), ROOT_URLCONF='obsdb.urls')

    import django
    django.setup()

    global MeasuredCCD
    import models
    MeasuredCCD = models.MeasuredCCD

    return settings

MeasuredCCD = None
