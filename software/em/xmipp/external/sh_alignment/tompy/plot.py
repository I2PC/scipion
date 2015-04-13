#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy

import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib as mlt
import matplotlib.patches as mpatches
from pylab import *

class TomImage:
    def __init__(
        self,
        data,
        x=-1,
        y=-1,
        z=-1,
        ):
        self.data = data

        self.x = x
        self.y = y
        self.z = z
        self.shape = np.array(data.shape) - 1

    def pltnormalize(self, p):
        # Normalize

        s = 3 * np.std(p)
        p[p > s] = s
        p[p < -s] = -s
        p = p - np.min(p)
        p = p / np.max(p)
        return p

    def getImage(self):
        if self.x != -1:
            p = self.data[self.x, :, :]
        elif self.y != -1:
            p = self.data[:, self.y, :]
        elif self.z != -1:
            p = self.data[:, :, self.z]
        else:
            raise Exception('Must provide a valide slice number!')


        return self.pltnormalize(p)

        # return self.plnormalize(np.swapaxes(p,1,0))

    def inc(self):#Incease the slice number
        if self.x > -1 and self.x < self.shape[0]:
            self.x += 1
        elif self.y > -1 and self.y < self.shape[1]:
            self.y += 1
        elif self.z > -1 and self.z < self.shape[2]:
            self.z += 1
        else:
            pass

    def dec(self):#Decrease the slice number
        if self.x > 0:
            self.x -= 1
        elif self.y > 0:
            self.y -= 1
        elif self.z > 0:
            self.z -= 1
        else:
            pass

    def getcaption(self):
        if self.x > -1:
            temp = 'The slice is X=' + str(self.x)
        elif self.y > -1:
            temp = 'The slice is Y=' + str(self.y)
        elif self.z > -1:
            temp = 'The slice is Z=' + str(self.z)
        else:
            pass
        return temp

def show2D(
    data,
    x=-1,
    y=-1,
    z=-1,
    filename = None
    ):
    """Display the data of certain slice.
    'up' and 'right' key on keyboard will increse the slice number, 'down' and 'left' will decrease the slice number.

    @param data: data to be shown.
    @param x: slice nr. of x
    @param y: slice nr. of y
    @param z: slice nr. of z
    """
    if len(data.shape) != 3:
        raise Exception("Data should be 3D!")

    # if nothing is specified use the center slice along z
    if x == -1 and y == -1 and z == -1:
        z = data.shape[2]/2
    
    img = TomImage(data, x, y, z)
    p = img.getImage()
    
    import matplotlib.cm as cm
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    import matplotlib as mlt

    fig = plt.figure()
    
    def on_key(event):
        if event.key == 'up' or event.key == 'right':
            img.inc()
        elif event.key == 'down' or event.key == 'left':
            img.dec()
        else:
            pass

        p = img.getImage()
        txt = img.getcaption()
        plt.clf()
        fig.text(0.5, 1, txt, ha='center', va='top')
        temp = plt.imshow(p, cmap=cm.gray, origin='upper')

        fig.canvas.draw()

    #
    cid = fig.canvas.mpl_connect('key_release_event', on_key)
    txt = img.getcaption()
    fig.text(0.5, 1, txt, ha='center', va='top')
    plt.imshow(p, cmap=cm.gray, origin='upper')

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, bbox_inches=0)

