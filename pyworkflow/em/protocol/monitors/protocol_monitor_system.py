# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.protocol.params as params
from protocol_monitor import ProtMonitor, Monitor
import sqlite3 as lite
import psutil
import time, sys

from pyworkflow import VERSION_1_1
from pyworkflow.gui.plotter import plt
from pyworkflow.protocol.constants import STATUS_RUNNING, STATUS_FINISHED
from pyworkflow.protocol import getProtocolFromDb
from pyworkflow.em.plotter import EmPlotter


#cuda related stuff
import platform
import ctypes

class DeviceProp(ctypes.Structure):
    _fields_ = [
         ("name", 256*ctypes.c_char), #  < ASCII string identifying device
         ("totalGlobalMem", ctypes.c_size_t), #  < Global memory available on device in bytes
         ("sharedMemPerBlock", ctypes.c_size_t), #  < Shared memory available per block in bytes
         ("regsPerBlock", ctypes.c_int), #  < 32-bit registers available per block
         ("warpSize", ctypes.c_int), #  < Warp size in threads
         ("memPitch", ctypes.c_size_t), #  < Maximum pitch in bytes allowed by memory copies
         ("maxThreadsPerBlock", ctypes.c_int), #  < Maximum number of threads per block
         ("maxThreadsDim", 3*ctypes.c_int), #  < Maximum size of each dimension of a block
         ("maxGridSize", 3*ctypes.c_int), #  < Maximum size of each dimension of a grid
         ("clockRate", ctypes.c_int), #  < Clock frequency in kilohertz
         ("totalConstMem", ctypes.c_size_t), #  < Constant memory available on device in bytes
         ("major", ctypes.c_int), #  < Major compute capability
         ("minor", ctypes.c_int), #  < Minor compute capability
         ("textureAlignment", ctypes.c_size_t), #  < Alignment requirement for textures
         ("deviceOverlap", ctypes.c_int), #  < Device can concurrently copy memory and execute a kernel
         ("multiProcessorCount", ctypes.c_int), #  < Number of multiprocessors on device
         ("kernelExecTimeoutEnabled", ctypes.c_int), #  < Specified whether there is a run time limit on kernels
         ("integrated", ctypes.c_int), #  < Device is integrated as opposed to discrete
         ("canMapHostMemory", ctypes.c_int), #  < Device can map host memory with cudaHostAlloc/cudaHostGetDevicePointer
         ("computeMode", ctypes.c_int), #  < Compute mode (See ::cudaComputeMode)
         ("__cudaReserved", 36*ctypes.c_int)]


    def __str__(self):
        return """NVidia GPU Specifications:
    Name: %s
    Total global mem: %f Gb
    Shared mem per block: %i
    Registers per block: %i
    Warp size: %i
    Mem pitch: %i
    Max threads per block: %i
    Max treads dim: (%i, %i, %i)
    Max grid size: (%i, %i, %i)
    Total const mem: %i
    Compute capability: %i.%i
    Clock Rate (GHz): %f
    Texture alignment: %i
""" % (self.name, self.totalGlobalMem/(1073741824.), self.sharedMemPerBlock,
       self.regsPerBlock, self.warpSize, self.memPitch,
       self.maxThreadsPerBlock,
       self.maxThreadsDim[0], self.maxThreadsDim[1], self.maxThreadsDim[2],
       self.maxGridSize[0], self.maxGridSize[1], self.maxGridSize[2],
       self.totalConstMem, self.major, self.minor,
       float(self.clockRate)/1.0e6, self.textureAlignment)

class Cuda(object):
    #instead of using pyCuda I am going to access C library through ctypes
    #may be this should be change in the future
    cudaSuccess = 0
    errorDict = {
        1: 'MissingConfigurationError',
        2: 'MemoryAllocationError',
        3: 'InitializationError',
        4: 'LaunchFailureError',
        5: 'PriorLaunchFailureError',
        6: 'LaunchTimeoutError',
        7: 'LaunchOutOfResourcesError',
        8: 'InvalidDeviceFunctionError',
        9: 'InvalidConfigurationError',
        10: 'InvalidDeviceError',
        11: 'InvalidValueError',
        12: 'InvalidPitchValueError',
        13: 'InvalidSymbolError',
        14: 'MapBufferObjectFailedError',
        15: 'UnmapBufferObjectFailedError',
        16: 'InvalidHostPointerError',
        17: 'InvalidDevicePointerError',
        18: 'InvalidTextureError',
        19: 'InvalidTextureBindingError',
        20: 'InvalidChannelDescriptorError',
        21: 'InvalidMemcpyDirectionError',
        22: 'AddressOfConstantError',
        23: 'TextureFetchFailedError',
        24: 'TextureNotBoundError',
        25: 'SynchronizationError',
        26: 'InvalidFilterSettingError',
        27: 'InvalidNormSettingError',
        28: 'MixedDeviceExecutionError',
        29: 'CudartUnloadingError',
        30: 'UnknownError',
        31: 'NotYetImplementedError',
        32: 'MemoryValueTooLargeError',
        33: 'InvalidResourceHandleError',
        34: 'NotReadyError',
        0x7f: 'StartupFailureError',
        10000: 'ApiFailureBaseError'}

    def __init__(self):
        self.getCudaLib()

    #if there is a cuda library _libcudart will NOT be None
    def getCudaLib(self):
        sys.stderr.write("getCudaLib_0\n")
        try:
            if platform.system() == "Microsoft":
                self._libcudart = ctypes.windll.LoadLibrary('cudart.dll')
            elif platform.system()=="Darwin":
                self._libcudart = ctypes.cdll.LoadLibrary('libcudart.dylib')
            else:
                sys.stderr.write("getCudaLib_1\n")
                self._libcudart = ctypes.cdll.LoadLibrary('libcudart.so')
                sys.stderr.write("getCudaLib_2: " + str(self._libcudart) + "\n")
            self._libcudart_error = None
        except OSError, e:
            self._libcudart_error = e
            self._libcudart = None
        sys.stderr.write("getCudaLib_3: " + str(self._libcudart)+"\n")
        exit()

    def getDriverVersion(self):
        if self._libcudart is None: return  None
        version = ctypes.c_int()
        self._libcudart.cudaDriverGetVersion(ctypes.byref(version))
        v = "%d.%d" % (version.value//1000,
                       version.value%100)
        return v

    def getRuntimeVersion(self):
        if self._libcudart is None: return  None
        version = ctypes.c_int()
        self._libcudart.cudaRuntimeGetVersion(ctypes.byref(version))
        v = "%d.%d" % (version.value//1000,
                       version.value%100)
        return v

    def _checkCudaStatus(self, status):
        if status != self.cudaSuccess:
            eClassString = self.errorDict[status]
            # Get the class by name from the top level of this module
            eClass = globals()[eClassString]
            raise eClass()

    def cudaGetDeviceCount(self):
        if self._libcudart is None: return  0
        deviceCount = ctypes.c_int()
        status = self._libcudart.cudaGetDeviceCount(ctypes.byref(deviceCount))
        self._checkCudaStatus(status)
        return deviceCount.value

    def _checkDeviceNumber(self, device):
        assert isinstance(device, int), "device number must be an int"
        assert device >= 0, "device number must be greater than 0"
        assert device < 2**8-1, "device number must be < 255"

    def getDeviceProperties(self, device):
        if self._libcudart is None: return  None
        self._checkDeviceNumber(device)
        props = DeviceProp()
        status = self._libcudart.cudaGetDeviceProperties(ctypes.byref(props), device)
        self._checkCudaStatus(status)
        return props

#end cuda stuff

SYSTEM_LOG_SQLITE = 'system_log.sqlite'


class ProtMonitorSystem(ProtMonitor):
    """ check CPU, mem and IO usage.
    """
    _label = 'system_monitor'
    _version = VERSION_1_1

    def __init__(self, **kwargs):
        ProtMonitor.__init__(self, **kwargs)
        self.dataBase = 'log.sqlite'
        self.tableName = 'log'

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):    
        ProtMonitor._defineParams(self, form)
        form.addParam('cpuAlert', params.FloatParam,default=101,
              label="Raise Alarm if CPU > XX%",
              help="Raise alarm if memory allocated is greater "
                   "than given percentage")
        form.addParam('memAlert', params.FloatParam,default=101,
              label="Raise Alarm if Memory > XX%",
              help="Raise alarm if cpu allocated is greater "
                   "than given percentage")
        form.addParam('swapAlert', params.FloatParam,default=101,
              label="Raise Alarm if Swap > XX%",
              help="Raise alarm if swap allocated is greater "
                   "than given percentage")

        form.addParam('monitorTime', params.FloatParam,default=300,
              label="Total Logging time (min)",
              help="Log during this interval")

        ProtMonitor._sendMailParams(self, form)
        group = form.addGroup('GPU')
        group.addParam('doGpu', params.BooleanParam, default=False,
                       label="Check GPU",
                       help="Set to true if you want to monitor the GPU")
        group.addParam('gpusToUse', params.StringParam, default='0',
                          label='Which GPUs to use:', condition='doGpu',
                          help='Providing a list of which GPUs '
                               '(0,1,2,3, etc). Default is monitor GPU 0 only')


    #--------------------------- STEPS functions ------------------------------

    def monitorStep(self):
        self.createMonitor().loop()

    def createMonitor(self):
        protocols = []
        for protPointer in self.inputProtocols:
            prot = protPointer.get()
            prot.setProject(self.getProject())
            protocols.append(prot)

        sysMonitor = MonitorSystem(protocols,
                                   workingDir=self.workingDir.get(),
                                   samplingInterval=self.samplingInterval.get(),
                                   monitorTime=self.monitorTime.get(),
                                   email=self.createEmailNotifier(),
                                   stdout=True,
                                   cpuAlert=self.cpuAlert.get(),
                                   memAlert=self.memAlert.get(),
                                   swapAlert=self.swapAlert.get(),
                                   doGpu=self.doGpu.get(),
                                   gpusToUse = self.gpusToUse.get())
        return sysMonitor

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        #TODO if less than 20 sec complain
        return []  # no errors

    def _summary(self):
        return ['Stores CPU, memory and swap ussage in percentage']

    def _methods(self):
        return []


class MonitorSystem(Monitor):
    """ This will will be monitoring a CTF estimation protocol.
    It will internally handle a database to store produced
    CTF values.
    """
    def __init__(self, protocols, **kwargs):
        Monitor.__init__(self, **kwargs)
        self.protocols = protocols
        self.cpuAlert = kwargs['cpuAlert']
        self.memAlert = kwargs['memAlert']
        self.swapAlert = kwargs['swapAlert']
        self._dataBase = kwargs.get('dbName', SYSTEM_LOG_SQLITE)
        self._tableName = kwargs.get('tableName', 'log')
        self.doGpu = kwargs['doGpu']
        if self.doGpu:
            #get Gpus to monitor
            self.gpusToUse = [int(n) for n in (kwargs['gpusToUse']).split()]
            #init GPU
            cuda = Cuda()
            #print some general info
            #yhis is helpfull for debuging but I do not think should be here

            sys.stderr.write("Driver version: %s\n" % cuda.getDriverVersion())
            sys.stderr.write("Runtime version: %s\n" % cuda.getRuntimeVersion())
            nn = cuda.cudaGetDeviceCount()
            sys.stderr.write("Device count: %s\n" % nn)
            for ii in self.gpusToUse:
                props = cuda.getDeviceProperties(ii)
                sys.stderr.write("\nDevice %d:\n" % ii)
                #sys.stderr.write(props.__str__())
                #sys.stderr.write(props)
                #for f_name, f_type in props._fields_:
                #    attr = props.__getattribute__(f_name)
                #    sys.stderr.write( "  %s: %s\n" % (f_name, attr))
            sys.sderr.write("Finished printing devices")
        else:
            self.gpusToUse = None

        self.conn = lite.connect(os.path.join(self.workingDir, self._dataBase),
                                 isolation_level=None)
        self.cur = self.conn.cursor()

    def _getUpdatedProtocol(self, prot):
        prot2 = getProtocolFromDb(prot.getProject().path,
                                  prot.getDbPath(),
                                  prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def warning(self, msg):
        self.notify("Scipion System Monitor WARNING", msg)

    def initLoop(self):
        self._createTable()
        psutil.cpu_percent(True)
        psutil.virtual_memory()

    def step(self):
        cpu = psutil.cpu_percent(interval=0)
        mem = psutil.virtual_memory()
        swap = psutil.swap_memory()

        if self.cpuAlert < 100 and cpu > self.cpuAlert:
            self.warning("CPU allocation =%f." % cpu.percent)
            self.cpuAlert = cpu

        if self.memAlert < 100 and mem.percent > self.memAlert:
            self.warning("Memory allocation =%f." % mem.percent)
            self.memAlert = mem.percent

        if self.swapAlert < 100 and swap.percent > self.swapAlert:
            self.warning("SWAP allocation =%f." % swap.percent)
            self.swapAlert = swap.percent

        sql = """INSERT INTO %s(mem,cpu,swap) VALUES(%f,%f,%f);""" % (
            self._tableName, mem.percent, cpu, swap.percent)

        try:
            self.cur.execute(sql)
        except Exception as e:
            print("ERROR: saving one data point (monitor). I continue")

        # Return finished = True if all protocols have finished
        return all(self._getUpdatedProtocol(prot).getStatus() != STATUS_RUNNING
                   for prot in self.protocols)

    def _createTable(self):
        self.cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                                id INTEGER PRIMARY KEY AUTOINCREMENT,
                                timestamp DATE DEFAULT (datetime('now','localtime')),
                                cpu FLOAT,
                                mem FLOAT,
                                swap FLOAT)
                                """ % self._tableName)

    def getData(self):
        cur = self.cur
        # I guess this can be done in a single call
        # I am storing the first meassurement
        cur.execute("select julianday(timestamp)  from %s where id=1" %
                    self._tableName )
        initTime = cur.fetchone()[0]
        cur.execute("select timestamp  from %s where id=1" % self._tableName)
        initTimeTitle = cur.fetchone()[0]
        cur.execute("select (julianday(timestamp) - %f)*24  from %s" %
                    (initTime,self._tableName) )
        idValues = [r[0] for r in cur.fetchall()]

        def get(name):
            try:
                cur.execute("select %s from %s" % (name, self._tableName))
            except Exception as e:
                print("ERROR readind data (plotter). I continue")
            return [r[0] for r in cur.fetchall()]

        data = {'initTime': initTime,
                'initTimeTitle': initTimeTitle,
                'idValues': idValues,
                'cpu': get('cpu'),
                'mem': get('mem'),
                'swap': get('swap'),
                }
        #conn.close()
        return data


from pyworkflow.viewer import ( DESKTOP_TKINTER, WEB_DJANGO, Viewer)
from matplotlib import animation


class ProtMonitorSystemViewer(Viewer):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'system monitor'
    _targets = [ProtMonitorSystem]

    def _visualize(self, obj, **kwargs):
        return [SystemMonitorPlotter(obj.createMonitor())]


class SystemMonitorPlotter(EmPlotter):
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'System Monitor'
    _targets = [ProtMonitorSystem]

    def __init__(self, monitor):
        EmPlotter.__init__(self, windowTitle="system Monitor")
        self.monitor = monitor
        self.y2 = 0.
        self.y1 = 100.
        self.win = 250 # number of samples to be ploted
        self.step = 50 # self.win  will be modified in steps of this size
        self.createSubPlot(self._getTitle(),
                           "time (hours)", "percentage (or Mb for IO)")
        self.fig = self.getFigure()
        self.ax = self.getLastSubPlot()
        self.ax.margins(0.05)
        self.ax.grid(True)
        self.oldWin = self.win

        self.lines = {}
        self.init = True
        self.stop = False

    def _getTitle(self):
        return ("Use scrool wheel to change view window (win=%d)\n "
                "S stops, C continues plotting" % self.win)

    def onscroll(self, event):

        if event.button == 'up':
            self.win += self.step
        else:
            self.win -= self.step
            if self.win < self.step:
               self.win = self.step

        if self.oldWin != self.win:
            self.ax.set_title(self._getTitle())
            self.oldWin= self.win
        self.animate()
        EmPlotter.show(self)

    def press(self,event):
        sys.stdout.flush()
        if event.key == 'S':
            self.stop = True
            self.ax.set_title('Plot is Stopped. Press C to continue plotting')
        elif event.key == 'C':
            self.ax.set_title(self._getTitle())
            self.stop = False
            self.animate()
        EmPlotter.show(self)

    def has_been_closed(self,ax):
        fig = ax.figure.canvas.manager
        active_fig_managers = plt._pylab_helpers.Gcf.figs.values()
        return fig not in active_fig_managers

    def animate(self,i=0): #do NOT remove i
                                         #FuncAnimation add it as argument

        if self.stop:
            return

        data = self.monitor.getData()
        self.x = data['idValues']
        for k,v in self.lines.iteritems():
            self.y = data[k]

            lenght = len(self.x)
            imin = max(0,len(self.x) - self.win)
            xdata = self.x[imin:lenght]
            ydata = self.y[imin:lenght]
            v.set_data(xdata,ydata)

        self.ax.relim()
        self.ax.autoscale()
        self.ax.grid(True)
        self.ax.legend(loc=2).get_frame().set_alpha(0.5)

    def paint(self,labels):
        for label in labels:
            self.lines[label],=self.ax.plot([], [], '-',label=label)

        anim = animation.FuncAnimation(self.fig, self.animate,
                                       interval=self.monitor.samplingInterval * 1000)#miliseconds

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        EmPlotter.show(self)

    def show(self):
        self.paint(['mem','cpu', 'swap'])

