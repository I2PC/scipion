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
import pyworkflow.protocol.params as params
from pyworkflow.em import Protocol
from pyworkflow.utils import getFloatListFromValues
import numpy as np
from scipy import interpolate, signal
#import cv2  ### get center
#import scipy.io as sio ### cart2pol
import math

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
                      pointerClass='SetOfParticles',
                      label="Siemens star pattern",  
                      help="Siemens star pattern is the input image or the "
                           "stack of input images with different ZP "
                           "(and same angle) for calculating the MTF.\n"
                           "Note: There is a different MTF for each different "
                           "ZP position.")
        form.addParam('pixelSize', params.FloatParam,
                      label="Resolution (nm)",  
                      help='Resolution, or pixel size. e.g: 13.11')
        form.addParam('refNumber', params.IntParam, default=-1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Number of times the randomization is performed",
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
                      
        form.addParallelSection(threads=0, mpi=0)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        
        
        # supposed that SS file is setOfParticles to get Dimensions and etc.
        # all below code are a very first draft and need to be organized properly KHOSOOSAN dar return ha va Functions
        imSS = self.inputSiemensStar.get()
        dx = self.pixelSize
        nRef = self.refNumber
        orders = getFloatListFromValues(self.fractionOrders.get())
        ringPos = getFloatListFromValues(self.ringPosition.get())
        
        self._insertFunctionStep('getMTFfromSiemensStar', imSS, dx, nRef, 
                                orders, ringPos)
                 
        #self._insertFunctionStep('gatherResultsStep', ....)        
    #--------------------------- STEPS functions --------------------------------------------
    
    def getMTFfromSiemensStar(self, imSS, dx, nRef, orders, ringPos):        
        nOrders = len(orders)
        ringRange = np.zeros(len(ringPos)*2)
        for k in range(len(ringPos)):
            ringRange[2*k] = ringPos[k] - 250;
            ringRange[2*k+1] = ringPos[k] + 250

        ############ Defining dimensions ##########   #### az inja
        #iDimRaw = np.shape(imSS)
        #if len(iDimRaw) == 2:
        #    Nim = 1
        #    iDim = iDimRaw
        #elif len(iDimRaw) == 3:
            # In hdf5 the number of images is in the first index of the dims.
            # In hdf5, rows and cols are in the second and third index of dims. 
        #    Nim = iDimRaw[0];
        #    iDim = (iDimRaw[1], iDimRaw[2])
        #else:
        #    raise Exception('Wrong number of dimensions for input image imSS')    

        ##############################################   ### ta inja be validate montaghel shavad
        
        Nim = imSS.getSize()
        iDim = (imSS.getDim()[0], imSS.getDim()[1])
        
        
        
        # By default the reference image will be the central image in case 
        # that no input argument is given for nRef. 
        if nRef == -1:
            nRef = np.floor(Nim/2)
        ###########################################    
        
        
        
        
        ############ binarized SS FFT: taken as system input FFT ########
        #if Nim == 1 and len(iDimRaw) == 2:
        if Nim == 1:
            imTemp = imSS
        else:
            #imTemp = imSS[nRef,:,:]
            imTemp = imSS[nRef]

        # Define necessary objects to use methods from other classes
        #GetCenterObj = GetCenterSSClass() ####  def getCenter(self, im):   ######################################## estefadeye sahih barresi shavad
        #convert = ConvertCartPolar()  ######################################## estefadeye sahih barresi shavad

        #[ygc, xgc] = GetCenterObj.getCenter(imTemp)
        [ygc, xgc] = self._getCenter(imTemp)
        x1 = 0 - xgc
        x2 = iDim[1] - xgc
        y1 = 0 - ygc
        y2 = iDim[0] - ygc
        
        x = np.arange(x1, x2)
        y = np.arange(y1, y2)

        # Minimal radius using geo center.
        # We take the minimal radius that can be seen, with the complete ring.
        boundaries_mesh = [x1, x2, y1, y2]
        Rmax = round(min(np.absolute(boundaries_mesh))) - 16

        # Angular vector coordinates
        Ntheta = round(2*np.pi*(Rmax))
        thetaV = np.linspace(0, 2*np.pi, Ntheta)
        thetaV = thetaV[0:-1]
        Ntheta = len(thetaV)
       
        # Polar coordinates grid
        rv = np.arange(Rmax)
        [theta, ro] = np.meshgrid(thetaV, rv)

        rows = np.shape(theta)[0]
        cols = np.shape(theta)[1]   
        xxp = np.zeros((rows, cols))
        yyp = np.zeros((rows, cols))
        #(xxp, yyp) = convert.pol2cart(theta, ro)
        (xxp, yyp) = self._pol2cart(theta, ro)

        # Spatial resolution (Pixel size) at each radius
        dlr = (np.float32(rv)+1.0)/Rmax*dx
        # Frequency resolution at each radius
        dflr = 1.0/(dlr*Ntheta)
        # Frequency zero index
        fxc = np.floor(Ntheta/2.0)

        # Set last value according to the visibility of the last ring
        # Position in pixels of the void rings
        rRangeN = np.round(ringRange/dx) 
        pF = rRangeN
        idV = []
        for k in range(1, len(ringPos)):
            if Rmax > pF[2*k]:
                pF_partial = np.arange(pF[2*k-1], pF[2*k]+1)
                idV = np.concatenate([idV, pF_partial])
            else:
                pF_partial = np.arange(pF[2*k-1], Rmax)
                idV = np.concatenate([idV, pF_partial])
                break

        # Radius used to get the binary values in normalization
        rowL = idV[-1]
   
        # Polar decomposition
        polarStarRef = np.zeros((Rmax, Ntheta))
        sp = interpolate.RectBivariateSpline(y, x, imTemp)
        for i in range(np.shape(xxp)[0]):
            for j in range(np.shape(xxp)[1]):
                polarStarRef[i, j] = sp([xxp[i,j]], [yyp[i,j]])

        polarStarRefNorm = self._getpolarStarNorm(polarStarRef, rowL) #### estefadeye sahih check shavad

        # We compute the DFT of each row by means of the FFT algorithm.
        # I.e. a FFT vector for each ro radius.
        ifftshift_polarref = np.fft.ifftshift(polarStarRefNorm)
        polarRNFT = np.fft.fftshift(np.fft.fft(ifftshift_polarref))

        pRNFTEnd = np.absolute(polarRNFT[rowL-20,:])
        ref0 = abs(polarRNFT[rowL-20, fxc])

        # Positions of first diffraction orders
        pRNorder1Pos = np.where((pRNFTEnd > 0.2*ref0) & 
                                (pRNFTEnd < 0.8*ref0))[0]

        # Spatial Frequency index of firt diffraction order at each radius
        peakPeriodRN = pRNorder1Pos[1] - fxc
        # Spatial Frequency of first diffraction order at each radius
        # Multiplication of the pixel step in frequency by the number of pixels
        # where the first order peak is located.
        fx = dflr*peakPeriodRN

        # Sign inversion maybe considered in the future
        # signInvertedBW = 0 

        # Frequency positions to lineary interpolate the MTF values in nm.
        # 1024 is the number of pixels of Mistral images.
        #xdim = 1024
        xdim = 512 
        # I'm hardcoding the 512 because mtf2psf takes too much time
        # with 1024, because of the extend1d2d method.
        # Finally in the beamline 1024 should be used because TXM 
        # images are 1024*1024. For the tests I could subsample the
        # images of mistral to 512.
        
        print "**********starting mtfRef************"
        
        
        maxFx = 0.03
        fxlin = np.linspace(0, maxFx, xdim)
        idV = idV.astype(int)
        mtfref = np.zeros((xdim, nOrders))
        for kk in range(nOrders):
            row_notVoid = idV
            col_harmonic = fxc + peakPeriodRN*orders[kk]
            
            xx_base = fx[row_notVoid]*orders[kk]
            yy_base = np.absolute(polarRNFT[row_notVoid, col_harmonic])

            if xx_base[0] < xx_base[-1]:
                xxref = xx_base
                yyref = yy_base
            else:
                # We have to reverse the x and the y vectors in case that the
                # are not in the increasing order.
                xxref = np.fliplr([xx_base])[0]
                yyref = np.fliplr([yy_base])[0]                

            # The function to be used for interpolation is deduced here
            f = interpolate.interp1d(xxref, yyref, bounds_error=False)
            mtfref[:, kk] = f(fxlin)

        # We get the rest of frequencies from the last order
        idFmax = np.asarray(np.where(fxlin > fx[idV[0]]*orders[-1])).flatten()
        idFmax = idFmax[0]
        # We repeat the last value of the last order for the last nan values.
        # Higher frequencies (smaller radius) attenuation coefficients.
        mtfref[idFmax:, nOrders-1] = mtfref[idFmax-1, nOrders-1]
        
        
        print "**********end mtfRef************"
        print "mtfref = \n"
        print mtfref

        ##################################################################
        ########## Itereting over images for system output FFT ###########
        
        print "**********starting mtfB************"
        
        mtfb = np.zeros((Nim, xdim, nOrders))
        # Selecting image if input is a single image
        if Nim == 1 and len(iDimRaw) == 2:
            imTemp = imSS
        for k_img in range(Nim):
            # Selecting image if input is an image stack
            if len(iDimRaw) > 2:
                imTemp = imSS[k_img, :, :]

            [ygc, xgc] = self._getCenter(imTemp)
            x1 = 0 - xgc
            x2 = iDim[1] - xgc
            y1 = 0 - ygc
            y2 = iDim[0] - ygc
            
            x = np.arange(x1, x2)
            y = np.arange(y1, y2)

            # Minimal radius using geometric center.
            # We take the minimal radius that can be seen, 
            # with the complete ring.
            boundaries_mesh = [x1, x2, y1, y2]
            Rmax = round(min(np.absolute(boundaries_mesh))) - 16

            # Polar coordinates grid
            rv = np.arange(Rmax)
            [theta, ro] = np.meshgrid(thetaV, rv)

            rows = np.shape(theta)[0]
            cols = np.shape(theta)[1]   
            xxp = np.zeros((rows, cols))
            yyp = np.zeros((rows, cols))
            (xxp, yyp) = convert.pol2cart(theta, ro)

            # Polar composition
            polarStar = np.zeros((Rmax, Ntheta))
            sp = interpolate.RectBivariateSpline(y, x, imTemp)
            for i in range(np.shape(xxp)[0]):
                for j in range(np.shape(xxp)[1]):
                    polarStar[i, j] = sp([xxp[i,j]], [yyp[i,j]])

            # We compute the DFT of each row by means of the FFT algorithm.
            # I.e. a FFT vector for each ro radius.
            ifftshift_polar = np.fft.ifftshift(polarStar)
            polarFT = np.fft.fftshift(np.fft.fft(ifftshift_polar))
            order1Pos = pRNorder1Pos

            # Sign inversion maybe considered in the future
            # signInverted = 0;

            dlr = (np.float32(rv)+1.0)/Rmax*dx
            # Frequency resolution at each radius
            dflr = 1.0/(dlr*Ntheta)
            # Frequency zero index
            fxc = np.floor(Ntheta/2.0)
            peakPeriod = order1Pos[1] - fxc
            fx = dflr*peakPeriod
            pF = rRangeN

            idV = []
            for k in range(1, len(ringPos)):
                if Rmax > pF[2*k]:
                    pF_partial = np.arange(pF[2*k-1], pF[2*k]+1)
                    idV = np.concatenate([idV, pF_partial])
                else:
                    pF_partial = np.arange(pF[2*k-1], Rmax)
                    idV = np.concatenate([idV, pF_partial])
                    break

            idV = idV.astype(int)
            for k_order in range(nOrders-1):
                row_notVoid = idV
                col_harmonic = fxc + peakPeriod*orders[k_order]
                
                xx_base = fx[row_notVoid]*orders[k_order]
                yy_base = np.absolute(polarFT[row_notVoid, col_harmonic])

                if xx_base[0] < xx_base[-1]:
                    xxb = xx_base
                    yyb = yy_base
                else:
                    # We have to reverse the x and the y vectors in case that the
                    # are not in the increasing order.
                    xxb = np.fliplr([xx_base])[0]
                    yyb = np.fliplr([yy_base])[0]                

                # The function to be used for interpolation is deduced here
                f = interpolate.interp1d(xxb, yyb, bounds_error=False)
                mtfb[k_img, :, k_order] = f(fxlin)

            # We get the rest of frequencies from the last order
            idFmax = (np.asarray(np.where(fx<maxFx/orders[-1])).flatten())[0]
            idVTmp_partial = np.arange(idFmax, idV[0]-1)
            idVTmp_partial = idVTmp_partial.astype(int)
            idVTmp = np.concatenate([idVTmp_partial, idV])

            col_last_harmonic = fxc + peakPeriod*orders[-1]
            xx_base = fx[idVTmp]*orders[-1]
            yy_base = np.absolute(polarFT[idVTmp, col_last_harmonic])

            if xx_base[0] < xx_base[-1]:
                xxb_end = xx_base
                yyb_end = yy_base
            else:
                # We have to reverse the x and the y vectors in case that the
                # are not in the increasing order.
                xxb_end = np.fliplr([xx_base])[0]
                yyb_end = np.fliplr([yy_base])[0]                

            # The function to be used for interpolation is deduced here
            f = interpolate.interp1d(xxb_end, yyb_end, bounds_error=False)
            mtfb[k_img, :, -1] = f(fxlin)
            
        print "**********end mtfb************"
        print "mtfb = \n"
        print mtfb
    
        ########################################################################
        ### Deducing MTF from the calculated FFTs of input (ref) and outputs ###
        ######################## Unifying curves ###############################
        
        print "**********starting mtf************"

        mtf = np.zeros((xdim, Nim))
        for k_img in range(Nim):
            mtfTmp = np.zeros((xdim, nOrders))
            nanPos = np.zeros((xdim, nOrders))
            for k_order in range(nOrders):
                mtfb_part = mtfb[k_img, :, k_order]
                mtfref_part = mtfref[:, k_order]
                mtfTmp[:, k_order] = mtfb_part/mtfref_part
                nanPos[:, k_order] = np.isnan(mtfTmp[:, k_order])

                for i in range(xdim):
                    if nanPos[i, k_order] == 1:
                        mtfTmp[i, k_order] = 0.0

            # The beginning is set to NaN just to be discarded
            nanNormT = nOrders - np.sum(nanPos, 1)
            noNanPos = np.asarray(np.where(nanNormT != 0)).flatten()
            mtf[:, k_img] = mtfTmp[:, 0]
            mtf[0:noNanPos[0], k_img] = np.nan

            # Calculation of overlapped range limits between consecutive orders 
            # to combine with progressive weights.
            for k_order in range(nOrders-1):
                nanspos_one = np.where(nanPos[:, k_order] == 0)
                e1 = (np.asarray(nanspos_one)).flatten(nanspos_one)[-1]
                nanspos_two = np.where(nanPos[:, k_order+1] == 0)
                i2 = (np.asarray(nanspos_two)).flatten()[0]

                if k_order == 2 and orders[-1] == 5:
                    if e1-i2 > 10:
                        i2 = e1-50
                mtf[i2:, k_img] = mtf[i2:, k_img] + mtfTmp[i2:, k_order+1]

                if  e1>i2:
                    # The curves for different orders overlaps. 
                    # We apply to unify curves weights
                    wg = np.linspace(0, 1, e1-i2+1)
                    a_mult = mtfTmp[i2:e1+1, k_order]
                    b_mult = mtfTmp[i2:e1+1, k_order+1]
                    mtf[i2:e1+1, k_img] = a_mult*(1-wg) + b_mult*wg
                else:      
                    # The curves for different orders does not overlap
                    x_for_interp = fxlin[noNanPos]
                    y_for_interp = mtf[noNanPos, k_img]
                    f = interpolate.interp1d(x_for_interp, y_for_interp)
                    mtf[e1:i2, k_img] = f(fxlin[e1:i2])
                    
                    
        print "**********end mtf************"
        print "mtf = \n"
        print mtf

        # Output dictionary containing the mtf, the mtfref, the mtfb and the 
        # fxlin (named fx in the dictionary).
        mtf_out_dict = {}
        mtf_out_dict['fx'] = fxlin
        mtf_out_dict['mtfref'] = mtfref
        mtf_out_dict['mtfb'] = mtfb
        mtf_out_dict['mtf'] = mtf
        
        print "**********dic is ready************"

        return mtf_out_dict
        
            
    
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
    def _getpolarStarNorm(self, polarStar, nrows_radius=-1):

        rows = np.shape(polarStar)[0]
        cols = np.shape(polarStar)[1]
        if nrows_radius == -1:
            nrows_radius = rows
            ref_row = rows-1
        else:
            ref_row = nrows_radius-1    

        xdim_theta = cols
        pMean = np.average(polarStar, 1)
        imBW = np.zeros((rows, cols))
        for i in range(rows):
            for j in range(cols):
                if polarStar[i,j] > pMean[i]:
                    imBW[i,j] = 1

        dpos = np.diff(imBW[ref_row, :])
        ndpos = np.asarray(np.where(dpos < 0)).flatten()
        pdpos = np.asarray(np.where(dpos > 0)).flatten()

        if pdpos[0] < ndpos[0]:
            pdpos = np.roll(pdpos, -1)     

        ss = 5
        polarStarFixBW = np.zeros((xdim_theta))

        for kk in range(len(ndpos)):
            if ndpos[kk] < pdpos[kk]:
                vTmp = np.arange(ndpos[kk]+1, pdpos[kk]+1)
            else:      
                vTmp1 = np.arange(ndpos[kk]+1, xdim_theta)
                vTmp2 = np.arange(0, pdpos[kk]+1)
                vTmp = np.concatenate([vTmp1, vTmp2])            

            # We look for the average of a portion of the ref row of the polar
            # star. The average of the columns referenced by the indexes in vTmp
            # The ref row is a row close to the last row 
            # (near the external circle).
            avg = 0
            count = 0
            vTmp_partial = vTmp[ss:-1-ss]
            for j in range(len(vTmp_partial)):
                count = count + 1
                avg = avg + np.average(polarStar[ref_row, vTmp_partial[j]])
            if count != 0:
                avg = avg/count
            else:
                avg = 0
            for i in range(len(vTmp)):
                polarStarFixBW[vTmp[i]] = avg

            if kk < len(ndpos)-1:
                vTmp = np.arange(pdpos[kk]+1, ndpos[kk+1]+1)
            else:
                if pdpos[kk] < ndpos[1]:
                    vTmp = np.arange(pdpos[kk]+1, ndpos[0])
                else:
                    vTmp1 = np.arange(pdpos[-1]+1, xdim_theta)
                    vTmp2 = np.arange(0, ndpos[0]+1)
                    vTmp = np.concatenate([vTmp1, vTmp2])

            # We look again for the average of a portion of the ref row of the 
            # polar star. The average of the columns referenced by the indexes 
            # in vTmp
            avg = 0
            count = 0
            vTmp_partial = vTmp[ss:-1-ss]
            for j in range(len(vTmp_partial)):
                count = count + 1
                avg = avg + np.average(polarStar[ref_row, vTmp_partial[j]])
            if count != 0:
                avg = avg/count
            else:
                avg = 0
            for i in range(len(vTmp)):
                polarStarFixBW[vTmp[i]] = avg

        higher_than_pMean = []
        lower_than_pMean = []
        for i in range(len(polarStarFixBW)):
            if polarStarFixBW[i] > pMean[ref_row]:
                higher_than_pMean.append(polarStarFixBW[i])
            elif polarStarFixBW[i] < pMean[ref_row]:
                lower_than_pMean.append(polarStarFixBW[i])
                 
        maxm = np.average(higher_than_pMean)
        minm = np.average(lower_than_pMean)

        # We normalize the image thanks to the found maxm and minm.
        imBW = imBW*(maxm-minm) + minm;
        return imBW
    
    
    def _getCenter(self, im):
        """ The center of the image is found by rotating the image
            by 180 degrees, and see if the center corresponds or not.
            The objective of this function is to find the center of 
            the SS pattern in order to be able to transform the image
            to polar coordinates.
        """
        cv = [0, 0]
        #cv[0] = np.ceil((np.shape(im)[0]+1)/2);
        #cv[1] = np.ceil((np.shape(im)[1]+1)/2);
        
        cv[0] = np.ceil((im.getDim()[0]+1)/2);
        cv[1] = np.ceil((im.getDim()[1]+1)/2);

        # We take a ROI in the central zone approx.
        # 150 should be used, but as it is too slow, it has been shown that for 
        # the image 80 that is the one that I'm using (81 in Matlab), cSize = 70 
        # works still well.
        # I will put 20 to go fast, but afterward I have to put 70 again. And
        # eventually, put 150.
        cSize = 70
        i0 = im[cv[0]-cSize:cv[0]+cSize, cv[1]-cSize:cv[1]+cSize]
        i0 = i0 - np.mean(i0)
        roti0 = np.rot90(i0, 2)

        # Compute the crosscorrelation between the image and the image rotated
        # by 180 degrees (in order to find the center of the pattern later on).
        
        
        #method_cc = cv2.TM_CCORR
       
        
        
        # This crosscorrelation it is not very fast. We should try to optimize 
        # it. Maybe by using cv2 and a small template.
        cc = signal.correlate2d(i0, roti0)

        # Position of the crosscorrelation maximum value
        cc_abs = np.absolute(cc)
        [y_peak, x_peak] = np.unravel_index(cc_abs.argmax(), cc_abs.shape)

        ccc = [0, 0]
        ccc[0] = np.ceil((np.shape(cc)[0] + 1) / 2)
        ccc[1] = np.ceil((np.shape(cc)[1] + 1) / 2)

        yc = cv[0]-(ccc[0] - y_peak)/2
        xc = cv[1]-(ccc[1] - x_peak)/2

        return [yc, xc]
    
    
    
    

    def _convertImCart2Polar(self, imIn, center=[-1,-1], outSize=-1):
        """
        Convert an image to polar coordinates, taking into account the 
        coordinates origin placed at center = [xcenter ycenter].

        By default, the resolution is kept along radial direction. Angular
        resolution is such as the spatial resolution for the larger radius
        matches the radial resolution.

        convertImCart2Polar calculates the image in polar coordinates 
        with outSize = [rSize aSize], being rSize and aSize the number of pixels 
        for radial and angular dimensions, respectively.
        Original coding by Joaquin Oton, National Center for Biotechnology, CSIC
        (joaquin.oton@csic.es)
        Translation to Python by Marc Rosanes, ALBA CELLS Synchrotron 
        (mrosanes@cells.es)
        Version 1.0
        """

        iDim = [0, 0]
        if len(np.shape(imIn)) == 3:
            [nIm, iDim[0], iDim[1]]= np.shape(imIn)
        elif len(np.shape(imIn)) <= 2:
            nIm = 1
            [iDim[0], iDim[1]]= np.shape(imIn)

        if center == [-1,-1]:
            center[0] = math.ceil((iDim[0]+1)/2);
            center[1] = math.ceil((iDim[1]+1)/2);

        # There are difference between usage of image indexes in Matlab and in 
        # Python. Pay attention and see that hereafter we mix index 0 and 1.
        # By doing so we try to follow the Kino convention.
        # If results are not Ok at the end, have a look at these indexes in 
        # Python.
        # Bring the center of the SS pattern to the real center and compute the
        # mesh.  
        x1 = 0-center[0]
        x2 = iDim[1]-center[0]
        y1 = 0-center[1]
        y2 = iDim[0]-center[1]
        
        x = np.arange(x1, x2)
        y = np.arange(y1, y2)

        # Minimal radius using geo center.
        # Why it is needed this radius? A possible answer is that we need the
        # minimal radius that can be seen, with the complete ring.
        boundaries_mesh = [x1, x2, y1, y2]
        Rmax = round(min(np.absolute(boundaries_mesh)))
        if outSize == -1:
            # Radius vector
            NR = Rmax + 1
            Ntheta = math.ceil(2*math.pi*Rmax)
        else:
            NR = outSize[0]
            Ntheta = outSize[1]

        # Radial coordinates
        rv = np.linspace(0, Rmax, NR)
        # Angular vector coordinates
        thetaV = np.linspace(0, 2*math.pi, Ntheta+1)
        thetaV = thetaV[0:-1]

        # Polar coordinates grid
        [theta, ro] = np.meshgrid(thetaV, rv)

        rows = np.shape(theta)[0]
        cols = np.shape(theta)[1]   
        xxp = np.zeros((rows, cols), dtype=np.float32)
        yyp = np.zeros((rows, cols), dtype=np.float32)
        xxp, yyp = self._pol2cart(theta, ro)

        # The interpolation part should be optimized by removing the for loops,
        # because it takes to long, but without the for's it goes out of 
        # memory. Try to optimize it in the future.
        # Interpolation
        if nIm == 1:
            polar_img = np.zeros((NR, Ntheta))
            sp = interpolate.RectBivariateSpline(y, x, imIn)
            for i in range(np.shape(xxp)[0]):
                for j in range(np.shape(xxp)[1]):
                    polar_img[i][j] = sp([xxp[i,j]], [yyp[i,j]])

        else:
            # TODO: This else part has not been tested. To be done.
            for k in range(nIm):
                polar_img = np.zeros((NR, Ntheta, nIm))
                sp = interpolate.RectBivariateSpline(y, x, imIn[:,:,k])
                for i in range(np.shape(xxp)[0]):
                    for j in range(np.shape(xxp)[1]):
                        polar_img[i][j][k] = sp([xxp[i,j]], [yyp[i,j]])

        return polar_img
    
    ################################################################################################# in do aval dar tabe e bala eslah shavand sepas tabe e bala dar insert eslah shavad....
    def _pol2cart(self, theta, ro):
        if hasattr(theta, "__len__"):
            x = ro*np.cos(theta)
            y = ro*np.sin(theta)
        else:
            x = ro*math.cos(theta)
            y = ro*math.sin(theta)
        return (x, y)

    def _cart2pol(self, x, y):
        if hasattr(x, "__len__"):
            theta = np.arctan2(y, x)
            rho = np.sqrt(np.add(x**2 + y**2))
        else:
            theta = math.atan2(y, x)
            rho = math.sqrt(x**2 + y**2)
        return (theta, rho)