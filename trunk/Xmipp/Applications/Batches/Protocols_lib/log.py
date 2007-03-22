class logpyl:
    """The logpyl class implements basic logging functionality for Python programs.
    The log file contains the events logged to a logpyl object.
    It also allows to send mail confirmation
    Based on Logpyl - basic logfile objects for Python by john cole
    
    Author:Roberto Marabini, March 2007
"""

    def __init__(self,path='Logs',
                      name='logfile',
                      debug='off',
                      my_smtp_local_server='mail.cnb.uam.es'
                      ):
        "Initialize or open logs as log objects are instantiated."
        import sys        
        import os
        import socket

        self.my_smtp_local_server = my_smtp_local_server
        self.events = []    # list of events written to this log
        self.debug = debug  # flag to activate/deactive debugging messages
        
        # Try the log file
        self.name = name
        self.path = path
        self.logfile = name + '.log'
        self.ldf = path + '/' + self.logfile
        lfn = os.path.isfile(self.ldf)
        if ( lfn ):
            if ( debug == 'on' ):
                print 'DEBUG: Log file',self.logfile,'exists.'
            lfn = open(self.ldf,'r')
            for line in lfn.readlines():
                self.events.append(line.strip())
            lfn.close()
        else:
            if ( debug == 'on' ):
                print 'DEBUG: Log file',self.logfile,'does not exist.'
        # append a line with user, machine and date
        myusername = str(os.environ.get('USERNAME'))
        myhost = str(socket.gethostname())
        mypwd = str(os.environ.get('PWD'))
        event = "\n\n===============\n"
        event += "NEW LOG SESSION\n" 
        event += "===============\n" 
        event += myusername + '@' 
        event += myhost + ':' 
        event += mypwd + ' (' 
        event += self.datetime() + ')\n'
        self.events.append(event)

            
    def add(self, eventclass="note", message="Your message here"):
        "Compose a log entry from the elements passed to add() and append it to the list of events."
        import time
        event = eventclass + ' ' + message + ' (' + self.datetime() + ')'
        if ( self.debug == "on" ):
            print 'DEBUG: Adding', event, 'to log', self.name
        self.modified = time.asctime(time.localtime(time.time()))            
        self.events.append(event)
        return

    def close(self):
        "Close the log by writing all log entries to the proper file."
        import sys
        import os
        "Write the current version of the log to a file and free the variables used by the log."
        if ( self.debug == 'on' ):            
            print "DEBUG: Closing log", self.name
        # If self.path does not exist, create the directory for the logfiles.
        if ( not os.path.exists(self.path ) ):
            if ( self.debug == 'on' ):
                print 'DEBUG: Directory ',self.path,' does not exist.  I am creating it now.'
            try:                
                os.makedirs(self.path)
                if ( self.debug == 'on' ):
                    print 'DEBUG: Created log file directory',self.path
            except OSERROR:
                print 'ERROR: Could not create log file directory',self.path                

        # Make sure that the log entries are written.
        lfn = open(self.ldf, 'w+')
        for event in self.events:
            lfn.write(event+'\n')
        lfn.close()


    def datetime(self):
        "Generate the date/time stamp used in our log entries"
        import time
        datestamp = time.asctime(time.localtime(time.time()))
        return datestamp

    def printlog(self):
        print '\nPrinting log', self.name
        for event in self.events:
            print event
        print '\n'   


    ###########################
    # Send mail
    ##########################
    def my_send_mail(self,my_from,my_to,my_subject,my_message):
         import smtplib
         msg = ("From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n" % (my_from, my_to, my_subject))
         line = msg + my_message
         try:
             myserver=smtplib.SMTP(self.my_smtp_local_server)
             myserver.sendmail(my_from, my_to, line)
         except smtplib.SMTPException, error:
             print error
             ThisIsTheEnd("no error message")


    ###########################
    ##send a mail (text), with attachments
    ###########################
    def my_send_mail_attach(self,fro, to, subject, text, files=[]):
        import smtplib
        server=self.my_smtp_local_server;
        assert type(to)==list
        assert type(files)==list

        msg = MIMEMultipart()
        msg['From'] = fro
        msg['To'] = COMMASPACE.join(to)
        msg['Date'] = formatdate(localtime=True)
        msg['Subject'] = subject

        msg.attach( MIMEText(text) )

        for file in files:
            part = MIMEBase('application', "octet-stream")
            part.set_payload( open(file,"rb").read() )
            Encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"'
                           % os.path.basename(file))
            msg.attach(part)

        smtp = smtplib.SMTP(server)
        smtp.sendmail(fro, to, msg.as_string() )
        smtp.close()
 
if __name__ == '__main__':

    # create a new log or open an existing log with debugging turned on
    # (disable debugging messages by passing 'off' as the third parm to logpyl())
    mylog = logpyl('Logs','testlog','on')
    # add a couple of events to the log
    mylog.add("spam","Spam is US$1.95 per can.")
    mylog.add("eggs","Eggs are US$0.89 per dozen. sdfsdf dsfsdfs sdfsdf sdfsdfds sdfsdf sdfsdf sdfsdf sfsdf sdfsdf sdfsdf sdfs sdf  cc")
    # print the log entries
    mylog.printlog()
    # close the log
    mylog.close()
    
