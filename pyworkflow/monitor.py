from __future__ import print_function
import time
instances = {}
elapsedTime = {}


def getInstancesCount():
    return instances


def instanceCreated(cls):
    instances[cls] = instances.get(cls, 0) + 1


def instanceDestroyed(cls):
    instances[cls] = instances.get(cls, 0) - 1


def printInstances():
    for k, v in instances.items():
        print("%s: %d", k, v)


def waitForDebugger(seconds=20):
    print ("Waiting for debugger %d seconds." % seconds)
    from pyworkflow.utils import printTraceBack
    printTraceBack()

    while seconds > 0:
        time.sleep(1)
        # Set seconds to 0:
        # Execute this in the debugger:  seconds = 0
        print(str(seconds) + " seconds left.")
        seconds -= 1

# From: https://www.huyng.com/posts/python-performance-analysis

class Timer(object):

    indentation = 0
    blackList = []
    whiteList = []
    minms = 1

    def __init__(self, name=None, verbose=True):
        self.verbose = verbose
        self.name = name
        Timer.indentation += 1

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        Timer.indentation -= 1

        toPrint = (len(Timer.whiteList) == 0)

        if any(self.name.startswith(s) for s in Timer.blackList):
            toPrint = False

        if any(self.name.startswith(s) for s in Timer.whiteList):
            toPrint = True

        if self.verbose and toPrint and self.msecs > Timer.minms:
            print('{0}{1:10.0f}\tms\t{2}'.format('\t' * Timer.indentation, self.msecs, self.name))

# To monitor memory leaks:
# 1.- Check you have installed pympler: scipion run pip install pympler
# 2.- Uncomment lines below
# 3.- Add a with statement like:
#     with monitor.MemoryMonitor():
#         <code to analyze>
# from pympler.tracker import SummaryTracker
#
#
# class MemoryMonitor(object):
#
#     def __init__(self):
#         self.tracker = None
#
#     def __enter__(self):
#         self.tracker = SummaryTracker()
#         self.tracker.print_diff()
#         return self
#
#     def __exit__(self, *args):
#         self.tracker.print_diff()
