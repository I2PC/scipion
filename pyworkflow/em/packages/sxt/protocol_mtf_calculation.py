# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Joaquin Oton   (joton@cnb.csic.es)
# *              Marc Rosanes   (mrosanes@cells.es)
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


#import os
#import glob
#from xmipp import *
import pyworkflow.protocol.params as params
from pyworkflow.em import Protocol
from pyworkflow.utils import getFloatListFromValues
import numpy as np
import pickle
import matplotlib.pyplot as plt
#from scipy.ndimage.interpolation import zoom
#from scipy import interpolate, signal
#import math
#import pickle
#from scipy.ndimage.interpolation import zoom



#from pyworkflow.em.data import Volume
#from pyworkflow.em.protocol import ProtReconstruct3D
#from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
#from pyworkflow.utils import getFloatListFromValues
#from pyworkflow.utils.path import cleanPattern, cleanPath, copyFile
#import xmipp
#from pyworkflow.object import Float, String
#from math import sqrt
#from plotter import XmippPlotter


class XmippProtMtfCalculation(Protocol):
    """    
    
    ?????? NEEDS TO BE CHECKED
    
    GETMTFFROMSIEMENSSTAR: MTF profiles Calculation.
    GETMTFFROMSIEMENSSTAR(imSS,dx) calculates the MTF profile from the
    imSS image of a Siemens star pattern for a pixel size dx in nm.

    If imSS is a stack of images, then the central slice is chosen to
    obtain the binarized reference image to normalize.

    The inputs:
    imSS: is the input image or the stack of input images with different
    ZP (and same angle) for calculating the MTF. There is a different MTF
    for each different ZP position.
    dx: Resolution, or pixel size. e.g: 13.11

    The output mtf_out_dict is a dictionary with the following fields:

    - mtfref: the reference Fourier coefs curves used to normalize at
    each Fourier order shown as column vectors x # of orders.

    - mtfb: the raw Fourier coefs curves calculated for each image in
    imSS at each Fourier order, for different radius of the ss.

    - mtf: The MTF profiles obtained after normalize mtfb using mtfref
    and combine the curves for all the Fourier ordes. This is the more 
    relevant output.

    - fx: column vector of spatial frecuency coordinates corresponding
    to the column vectors oft mtf, mtfb and mtfref. Units are nm^-1.

    GETMTFFROMSIEMENSSTAR(imSS, dx, nRef, orders, ringPos) specifies
    parameters different from default values:

    - nRef: image index of imSS used to calculate de reference image. By
    default, nRef = -1  automatically selects central image (floor(N/2)+1)
    Number of the reference image in the stack (specific zpZ)
    If not set, the reference image will be the image of the center of 
    the stack.

    - orders: diffraction orders used to calculate MTF profiles. By
    default, orders = [1 3].

    - ringPos: center radius in nm of the void rings of Siemens star
    pattern to be removed from the calculation. Default values are set 
    to match Xradia manufactured SS pattern: ringPos = [1.5 3 6 12]*1e3.
          
    """
    
    _label = 'calculating MTF from SS'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSiemensStar', params.PointerParam, 
                      pointerClass='SetOfImages',
                      label="Siemens star pattern",  
                      help="Siemens star pattern is the input image or the "
                           "stack of input images with different ZP "
                           "(and same angle) for calculating the MTF.\n"
                           "Note: There is a different MTF for each different "
                           "ZP position.")
        #form.addParam('pixelSize', params.FloatParam,
        #              label="Resolution (nm)",  
        #              help='Resolution, or pixel size. e.g: 13.11')
        form.addParam('refNumber', params.IntParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Reference number",
                      help="image index of siemens star pattern used to "
                           "calculate de reference image.\n"
                           "By default, nRef = -1  automatically selects central "
                           "image (floor(N/2)+1).")
        form.addParam('fractionOrders', params.NumericListParam,
                      default="1 3",
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Orders",
                      help="Diffraction orders used to calculate MTF profiles.") 
        form.addParam('ringPosition', params.NumericListParam,
                      default="1.5e3 3e3 6e3 12e3",
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Ring positions (nm)",
                      help="Center radius in nm of the void rings of "
                           "Siemens star pattern to be removed from the "
                           "calculation. Default values are set to match "
                           "Xradia manufactured SS pattern") 
                      
        form.addParallelSection(threads=1, mpi=2)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        
        
        # supposed that SS file is setOfParticles to get Dimensions and etc.
        # all below code are a very first draft and need to be organized properly KHOSOOSAN dar return ha va Functions
        imSS = self.inputSiemensStar.get()
        #dx = self.pixelSize.get()
        nRef = self.refNumber.get()
        orders = getFloatListFromValues(self.fractionOrders.get())
        ringPos = getFloatListFromValues(self.ringPosition.get())
        
        self._insertFunctionStep('getMTFfromSiemensStar', imSS, nRef, 
                                orders, ringPos)
        
                 
        #self._insertFunctionStep('gatherResultsStep', ....)        
    #--------------------------- STEPS functions --------------------------------------------
    
    def getMTFfromSiemensStar(self, imSS, nRef, orders, ringPos):        
        
        
        print "nRef = ", nRef
        print "orders = ", orders
        print "ringPos = ", ringPos
        print "imSS = ", imSS
        
        fhInputImage = self._getExtraPath('input_imSS.mrc')
        imSS.writeStack(fhInputImage)
        
        
        from xpytools.readTomo import accessData
        accessDataObj = accessData()
        tomo = accessDataObj.readTomo(fhInputImage)
        imgs_ss = tomo[79:81, :, :]
        single_imgss = accessDataObj.readImg(fhInputImage, 80)
        print 'single image dimensions=\n' , np.shape(single_imgss)
        print 'set of images dimensions=\n' , np.shape(tomo)
        print 'partial set of images dimensions=\n' , np.shape(imgs_ss)
        
        
        
        from xpytools.getResolutionfromSiemensStar import MTFgetResolution
        resolutionObj = MTFgetResolution()
        ss_img_resolution = resolutionObj.getResolutionfromSiemensStar(single_imgss, ringPos)
        print '#######   ******   #######'
        print "Resolution = \n" , ss_img_resolution
        
        
        
        
        dx = ss_img_resolution[0]
        print "dx=", dx ,"\n"
        from xpytools.getMTFfromSiemensStar import MTFfromSiemensStar
        MTFObj = MTFfromSiemensStar()
        #mtf_out = MTFObj.getMTFfromSiemensStar(single_imgss, dx, nRef, orders, ringPos) 
        #mtf_out = MTFObj.getMTFfromSiemensStar(tomo, dx, nRef, orders, ringPos)
        mtf_out = MTFObj.getMTFfromSiemensStar(imgs_ss, dx, nRef, orders, ringPos)
        pickle.dump(mtf_out, open("mtfoutputsingleimg_512.p", "wb")) #### 512 in the code should be checked and replace with the proper number....ask from Kino
        mtf_out = pickle.load(open("mtfoutputsingleimg_512.p", "rb"))
        fx = mtf_out['fx']
        mtfb = mtf_out['mtfb']
        mtfref = mtf_out['mtfref']
        mtf = mtf_out['mtf']
        print '#######   ******   #######'
        print "mtf dimension = " , np.shape(mtf)
        
        """
        plt.plot(fx, mtf[:,0],'-')
        plt.xlabel('Frequency (1/nm)')
        plt.ylabel('MTF')
        plt.title('MTF basd on a single image')
        plt.grid()
        plt.show()
        """
    
        
        from xpytools.mtf2psf import MTF2PSFClass
        mtf2psfObj = MTF2PSFClass()
        psfdict = mtf2psfObj.mtf2psf(mtf, fx, dx, fov=400, fc=-1) ## fov and fc should be an input value  or a fixed value like this in the main code??
        pickle.dump(psfdict, open("psf_xdim512_dx10_fov400_singleimg.p", "wb"))
        psfdict = pickle.load(open("psf_xdim512_dx10_fov400_singleimg.p", "rb"))
        psfarray = psfdict['psf']
        print '#######   ******   #######'
        print "psfarray dimension = " , np.shape(psfarray)
        
        """
        psfimg = psfarray[:,:,0]
        im = plt.imshow(psfimg, cmap="hot")
        plt.colorbar(im, orientation='horizontal')
        plt.show()
        """
            
    
    #def gatherResultsStep(self, ...):
         
    #--------------------------- INFO functions -------------------------------------------- 
    
#    def _summary(self):
#        """ Should be overriden in subclasses to 
#        return summary message for NORMAL EXECUTION. 
#        """
#              
#        msg = []
#        msg.append()        
#        msg.append()
#        return msg
#    
#    def _methods(self):
#        messages = []
#        messages.append('Joton')
#        return messages
#    
#    def _citations(self):
#        return ['Joton']
#    
#    def _validate(self):   ########################################## ba tavajoh be code haye Kino anjam shavad #####################################
#        errors=[]
#        if :
#            errors.append() 
#        if :
#            errors.append()
#        return errors 
             
    #--------------------------- UTILS functions --------------------------------------------
    