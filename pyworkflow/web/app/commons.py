import os

from pyworkflow.web.pages import settings

def staticPath(*paths):
    return os.path.join(settings.STATIC_URL, *paths)

