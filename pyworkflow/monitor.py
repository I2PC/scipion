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
        print "%s: %d", k, v


# From: https://www.huyng.com/posts/python-performance-analysis

class Timer(object):

    indentation = 0
    blackList = []
    whiteList = ['project.getRuns']
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
            print '{0} {1:10.0f} ms : {2}'.format(' ' * Timer.indentation, self.msecs, self.name)