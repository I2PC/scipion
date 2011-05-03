#!/usr/bin/env python

import smtplib
import config

def mail(toaddrs, fromaddr, subject, TEXT):
    import smtplib
    SERVER = config.SERVER
    #toaddrs  must be a list

    # Prepare actual message
    message = """\
From: %s
To: %s
Subject: %s

%s
    """ % (fromaddr, ", ".join(toaddrs), subject, TEXT)

    # Send the mail
    server = smtplib.SMTP(SERVER)
    server.sendmail(fromaddr, toaddrs, message)
    server.quit()

#message   = "my message"
#fromaddr = "xmipp@bioweb.cnb.csic.es"
#toaddrs  = "roberto.marabini.cnb@gmail.com"
#subject = "xmipp Compilation"
#mail(toaddrs, fromaddr, subject, message)

