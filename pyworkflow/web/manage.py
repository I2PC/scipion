#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":

    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "pages.settings")

    from django.core.management import execute_from_command_line


    # import time
    #
    # print ("Sleeping...")
    # time.sleep(10)
    # print ("Awake")

    execute_from_command_line(sys.argv)
