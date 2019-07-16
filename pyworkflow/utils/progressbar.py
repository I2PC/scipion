# **************************************************************************
# *
# * Authors:      Roberto Marabini (roberto@cnb.csic.es) [1]
# *               J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] SciLifeLab, Stockholm University
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

from __future__ import print_function
import sys
import re


class ProgressBar(object):
    """ Text progress bar class for Python.

    A text progress bar is typically used to display the progress of a long
    running operation, providing a visual cue that processing is underway.

    Example:
        N = 1000
        pb = ProgressBar(N, fmt=ProgressBar.FULL)
        pb.start()
        for x in xrange(N):
            pb.update(x+1)
            sleep(0.1)
        pb.finish()
    """
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'
    # scipion uses variable fonts so the bar size changes
    # since the space widt is different from the = width
    NOBAR = '%(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'
    OBJID = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) (objectId=%(objectId)d)'
    DOT = '.'

    def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',
                 output=None, extraArgs=None):
        """
        Create a new ProgressBar object.
        :param total: The total amount that will be running the progress bar.
            The value in the update() method can no go greater than this value.
        :param width: progress bar width (without the percentage and number of
            iterations loop)
        :param fmt: predefined format string, so far DEFAULT, FULL, OBJID and
            DOT are defined.
        :param symbol: progress bar is made with this symbol
        :param output:
        :param extraArgs: Additional arguments that can be passed to be used
            the fmt format. (e.g extraArgs={'objectId': 1} for fmt=OBJID
        """
        if len(symbol) != 1:
            raise Exception("Symbol should be only 1 character length. ")

        self._total = total
        self._width = width
        self._symbol = symbol
        self._output = output or sys.stdout
        self._current = -1
        self._extraArgs = extraArgs or {}
        self._fmt = fmt
        self._directPrint = fmt == self.DOT

        if not self._directPrint:
            # This line computes the number of digits
            # in total and rewrites the fmt string
            # so if total = 100, %d is converted to %3d
            self._fmt = re.sub(r'(?P<name>%\(.+?\))d',
                               r'\g<name>%dd' % len(str(total)), fmt)

    def __getStr(self):
        """ Internal function to return the current string value.
        It should be called after the value has being set.
        """
        if self._directPrint:  # print just a dot
            return self._fmt if self._current else ''

        percent = self._current / float(self._total)
        size = int(self._width * percent)
        remaining = self._total - self._current
        bar = '[' + self._symbol * size + ' ' * (self._width - size) + ']'

        args = {
            'total': self._total,
            'bar': bar,
            'current': self._current,
            'percent': percent * 100,
            'remaining': remaining,
        }
        args.update(self._extraArgs)

        return '\r' + self._fmt % args

    def start(self):
        """ Print empty progress bar. """
        self.update(0)

    def update(self, value):
        """
        Update the current value and print the progress.
        :param value: New value, should be greater than the previous
            value and less or equal the total value/
        :return:
        """
        if value < 0 or value <= self._current or value > self._total:
            raise Exception("Incorrect value provided. It should be greater "
                            "than previous value and between 0 and total. ")
        self._current = value
        self._output.write(self.__getStr())
        self._output.flush()

    def finish(self, printNewLine=True):
        """ Finalize the progress and
        print last update with 100% complete message """
        if self._current < self._total:
            self.update(self._total)
        # print a new line
        if printNewLine:
            self._output.write("\n")
