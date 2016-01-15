
def django_setup():
    import os
    from django.conf import settings
    basedir = os.path.dirname(os.path.dirname(__file__))
    settings.configure(INSTALLED_APPS=['obsdb'],
                       MIDDLEWARE_CLASSES=[],
                       DATABASES=dict(default=dict(
                           ENGINE='django.db.backends.sqlite3',
                           NAME=os.path.join(basedir,'obsdb','obsdb.sqlite3')))
        )

    import django
    django.setup()

    global MeasuredCCD
    import models
    MeasuredCCD = models.MeasuredCCD

    return settings

MeasuredCCD = None
