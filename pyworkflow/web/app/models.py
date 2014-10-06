from django.db import models

class Document(models.Model):
    docfile = models.FileField(upload_to='uploads/')
    