'''
Created on Sep 30, 2013

@author: adrian
'''

   

#machineName="glassfishprod.cnb.csic.es"
machineName="crunchy.cnb.csic.es"
usernameString="aquintana"
passwordString="q1w2e3r4"

import paramiko
import os
import select
import sys
import Xlib.support.connect as xlib_connect

def run(transport, session, command):
    def x11_handler(channel, (src_addr, src_port)):
        x11_fileno = channel.fileno()
        local_x11_channel = xlib_connect.get_socket(*local_x11_display[:3])
        local_x11_fileno = local_x11_channel.fileno()

        # Register both x11 and local_x11 channels
        channels[x11_fileno] = channel, local_x11_channel
        channels[local_x11_fileno] = local_x11_channel, channel

        poller.register(x11_fileno, select.POLLIN)
        poller.register(local_x11_fileno, select.POLLIN)

        transport._queue_incoming_channel(channel)

    def flush_out(channel):
        while channel.recv_ready():
            sys.stdout.write(channel.recv(4096))
        while channel.recv_stderr_ready():
            sys.stderr.write(channel.recv_stderr(4096))

    local_x11_display = xlib_connect.get_display(os.environ['DISPLAY'])

    channels = {}
    poller = select.poll()
    session_fileno = session.fileno()
    poller.register(session_fileno)

    session.request_x11(handler=x11_handler)
    session.exec_command(command)
    transport.accept()

    print "dondeestas"

    # event loop
    while not session.exit_status_ready():
        poll = poller.poll()
        if not poll: # this should not happen, as we don't have a timeout.
            break
        for fd, event in poll:
            if fd == session_fileno:
                flush_out(session)
            # data either on local/remote x11 channels/sockets
            if fd in channels.keys():
                sender, receiver = channels[fd]
                try:
                    receiver.sendall(sender.recv(4096))
                except:
                    sender.close()
                    receiver.close()
                    channels.remove(fd)

    flush_out(session)
    return session.recv_exit_status()

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh_client.connect('crunchy.cnb.csic.es', username='aquintana', password='q1w2e3r4')
    transport = ssh_client.get_transport()
    session = transport.open_session()
    run(transport, session, 'gedit')

#import os
#import select
#import sys
#
#import paramiko
#import Xlib.support.connect as xlib_connect
#
#import logging
#logging.basicConfig(level=logging.DEBUG)
#
#local_x11_display = xlib_connect.get_display(os.environ['DISPLAY'])
#local_x11_socket = xlib_connect.get_socket(*local_x11_display[:3])
#
#
#ssh_client = paramiko.SSHClient()
#ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#ssh_client.connect(machineName, username=usernameString, password=passwordString)
#transport = ssh_client.get_transport()
#session = transport.open_session()
#session.request_x11(single_connection=True)
#
#print "aki llegaste"
#session.exec_command('xp')
#print "yaki"
#x11_chan = transport.accept()
#print "ski no"
#
#session_fileno = session.fileno()
#x11_chan_fileno = x11_chan.fileno()
#local_x11_socket_fileno = local_x11_socket.fileno()
#
#poller = select.poll()
#poller.register(session_fileno, select.POLLIN)
#poller.register(x11_chan_fileno, select.POLLIN)
#poller.register(local_x11_socket, select.POLLIN)
#while not session.exit_status_ready():
#    poll = poller.poll()
#    if not poll: # this should not happen, as we don't have a timeout.
#        break
#    for fd, event in poll:
#        if fd == session_fileno:
#            while session.recv_ready():
#                sys.stdout.write(session.recv(4096))
#            while session.recv_stderr_ready():
#                sys.stderr.write(session.recv_stderr(4096))
#        if fd == x11_chan_fileno:
#            local_x11_socket.sendall(x11_chan.recv(4096))
#        if fd == local_x11_socket_fileno:
#            x11_chan.send(local_x11_socket.recv(4096))
#
#print 'Exit status:', session.recv_exit_status()
#while session.recv_ready():
#    sys.stdout.write(session.recv(4096))
#while session.recv_stderr_ready():
#    sys.stdout.write(session.recv_stderr(4096))
#session.close()
#    




print "machine", machineName