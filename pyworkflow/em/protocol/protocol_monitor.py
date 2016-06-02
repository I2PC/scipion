# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import sys

import pyworkflow.protocol.params as params
from protocol import EMProtocol



class ProtMonitor(EMProtocol):
    """ This is the base class for implementing 'Monitors', a special type
    of protocols intended to be used in the context of the Scipion-Box,
    where several protocols will be run in streaming and the monitor will
    be 'observing' their progress.
    """
    _label = 'monitor'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        form.addSection(label='Input')
        
        form.addParam('inputProtocols', params.MultiPointerParam,
                      label="Input protocols", important=True,
                      pointerClass='EMProtocol',
                      help="this protocol/s will be monitorized")

        form.addParam('samplingInterval', params.IntParam, default=60,
                      label="Sampling Interval (sec)",
                      help="Take one sample each SAmplinInteval seconds")

        #self._sendMailParams(form)

    def _sendMailParams(self, form):
        g = form.addGroup('Email settings')

        g.addParam('doMail', params.BooleanParam,
                      label="Enable Email notification?", default=False,
                      help="Allow monitors to notify via email.")

        g.addParam('emailFrom',params.StringParam, condition='doMail',
                   default = "from@from.fakeadress.com",
                   label='From',
                   help='Provide the sender address for notifications.')

        g.addParam('emailTo', params.StringParam, condition='doMail',
                   default = "to@to.fakeadress.com",
                   label='To',
                   help='Provide the destination address for notifications.')

        g.addParam('smtp', params.StringParam, condition='doMail',
                   default = "smtp.fakeadress.com",
                   label='SMTP Mail server',
                   help='Provide the address of SMTP mail server.')

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('monitorStep')

    #--------------------------- STEPS functions -------------------------------
    def monitorStep(self):
        pass

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        return []  # no errors

    def _summary(self):
        return []

    def _methods(self):
        return []

    def sendEMail(self, emailSubject, emailMessage):
        # Import smtplib for the actual sending function
        import smtplib
        # Import the email modules we'll need
        from email.mime.text import MIMEText

        msg = MIMEText(emailMessage)

        msg['Subject'] = emailSubject
        msg['From'] = self.emailFrom.get()
        msg['To'] = self.emailTo.get()

        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP(self.smtp.get())
        s.sendmail(self.emailFrom.get(), self.emailTo.get(), msg.as_string())
        s.quit()

    def createEmailNotifier(self):
        if getattr(self, 'doMail', False):
            email = EmailNotifier(self.smtp.get(),
                                  self.emailFrom.get(),
                                  self.emailTo.get())
        else:
            email = None

        return email


class Monitor():
    def __init__(self, **kwargs):
        # Where to store any data from this monitor
        self.workingDir = kwargs['workingDir']
        self.samplingInterval = kwargs.get('samplingInterval', None)
        self.monitorTime = kwargs.get('monitorTime', None)

        self._notifiers = []

        if kwargs.get('email', None) is not None:
            self._notifiers.append(kwargs['email'])

        if 'stdout' in kwargs:
            self._notifiers.append(PrintNotifier())

    def notify(self, title, message):
        for n in self._notifiers:
            n.notify(title, message)


class EmailNotifier():
    def __init__(self, smtpServer, emailFrom, emailTo):
        self._stmpServer = smtpServer
        self._emailFrom = emailFrom
        self._emailTo = emailTo

    def notify(self, title, message):
        # Import smtplib for the actual sending function
        import smtplib
        # Import the email modules we'll need
        from email.mime.text import MIMEText

        msg = MIMEText(message)

        msg['Subject'] = title
        msg['From'] = self._emailFrom
        msg['To'] = self._emailTo

        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP(self._smtpServer)
        s.sendmail(self._emailFrom, self._emailTo, msg.as_string())
        s.quit()


def PrintNotifier():
    def notify(title, message):
        print title, message
        sys.stdout.flush()