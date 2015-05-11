#!/usr/bin/env python

import smtplib
import mimetypes
import config
from email.mime.text import MIMEText
from email.MIMEMultipart import MIMEMultipart
from email.MIMEImage import MIMEImage
from email.MIMEBase import MIMEBase
from email.Encoders import encode_base64
from email import encoders

def mail(toaddrs, fromaddr, subject, TEXT, attach):
    import passwd
    #mail building
    if attach == None:
        msg = MIMEText(TEXT)
    else:
        msg = MIMEMultipart(TEXT)
        #Attaching the metadata
        for att in attach:
            _add_attachment(msg, att)
    msg['Subject'] = subject
    msg['From'] = fromaddr
    msg['To'] = toaddrs


    #authentication
    mailServer = smtplib.SMTP('smtp.gmail.com',587)
    mailServer.ehlo()
    mailServer.starttls()
    mailServer.ehlo()
    mailServer.login('biocompwebs@gmail.com', passwd.passwd)
    
    #sending
    mailServer.sendmail(fromaddr, toaddrs.split(","), msg.as_string())
    
    #close connection
    mailServer.close()
    
def _add_attachment(outer, filename):
    import sys
    import os
    ctype, encoding = mimetypes.guess_type(filename)
    if ctype is None or encoding is not None:
        # No guess could be made, or the file is encoded (compressed), so
        # use a generic bag-of-bits type.
        ctype = 'application/octet-stream'
    maintype, subtype = ctype.split('/', 1)
    fp = open(filename, 'rb')
    if maintype == 'text':
        # Note: we should handle calculating the charset
        msg = MIMEText(fp.read(), _subtype=subtype)
    elif maintype == 'image':
        msg = MIMEImage(fp.read(), _subtype=subtype)
    elif maintype == 'audio':
        msg = MIMEAudio(fp.read(), _subtype=subtype)
    else:
        msg = MIMEBase(maintype, subtype)
        msg.set_payload(fp.read())
        # Encode the payload using Base64
        encoders.encode_base64(msg)
    fp.close()
    # Set the filename parameter
    msg.add_header('Content-Disposition',
            'attachment',
            filename=os.path.basename(filename))
    outer.attach(msg)


