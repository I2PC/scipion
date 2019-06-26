from __future__ import print_function
import sys
import re

###############################################
# progressbar
# EXAMPLE:
#
# from time import sleep

# progress = ProgressBar(1000, fmt=ProgressBar.FULL)
# for x in xrange(progress.total):
#    progress.current += 1
#    print progress()
#    sleep(0.1)
# print progress.done()


class ProgressBar(object):
    """Text progress bar class for Python.

    A text progress bar is typically used to display the progress of a long
    running operation, providing a visual cue that processing is underway.
    """

    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'
    OBJID = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) (objectId=%(objectId)d)'
    DOT = '.'

    def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',
                 output=sys.stderr, objectId=None):
        # total = mamimum number of character to be written per line. Usually the current terminal width
        # width = progress bar width (without the percentange and number of iterations loop)
        # predefined format string, so far DEFAULT, FULL, OBJID and DOT are defined.
        # symbol: progress bar is made with this symbol

        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.objectId = objectId
        if fmt == self.DOT:
            self.directPrint = True
            self.fnt = self.DOT
        else:
            self.directPrint = False
            # This line computes the number of digits
            # in total and rewrites the fmt string
            # so if total = 100, %d is converted to %3d
            self.fmt = re.sub(r'(?P<name>%\(.+?\))d',
                r'\g<name>%dd' % len(str(total)), fmt)
        self.current = 0

    def __call__(self):
        # print just a dot
        if self.directPrint:
            return self.fnt
        else:  # print percentages
            percent = self.current / float(self.total)
            size = int(self.width * percent)
            remaining = self.total - self.current
            bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'

            args = {
                'total': self.total,
                'bar': bar,
                'current': self.current,
                'percent': percent * 100,
                'remaining': remaining,
                'objectId': self.objectId
            }
            return('\r' + self.fmt % args)

    def done(self):
        # print last update with 100% complete message
        self.current = self.total
        return(self())
