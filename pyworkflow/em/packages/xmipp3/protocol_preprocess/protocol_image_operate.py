# *****************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), Sep 2013
# * Ported to Scipion:
# *           Vahid Abrishami (vabrishami@cnb.csic.es), Oct 2014
# *
# * Refactored/Updated: Josue Gomez-Blanco (jgomez@cnb.csic.es), Jun 2016
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
# *****************************************************************************

from collections import OrderedDict

from pyworkflow.em.constants import ALIGN_NONE
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtOperateParticles, ProtOperateVolumes
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)

# Operands enum
OP_PLUS = 0
OP_MINUS = 1
OP_MULTIPLY = 2
OP_DIVIDE = 3
OP_MINIMUM = 4
OP_MAXIMUM = 5
OP_DOTPRODUCT = 6
OP_LOG = 7
OP_LOG10 = 8
OP_SQRT = 9
OP_ABS = 10
OP_POW = 11
OP_SLICE = 12
OP_COLUNM = 13
OP_ROW = 14
OP_RADIAL = 15
OP_RESET = 16

OP_CHOICES = OrderedDict() #[-1]*(OP_RESET+1)

OP_CHOICES[OP_PLUS]  = 'plus'
OP_CHOICES[OP_MINUS] = 'minus'
OP_CHOICES[OP_MULTIPLY] = 'multiply'
OP_CHOICES[OP_DIVIDE] = 'divide'
OP_CHOICES[OP_MINIMUM] = 'minimum'
OP_CHOICES[OP_MAXIMUM] = 'maximum'
OP_CHOICES[OP_DOTPRODUCT]  = 'dot product'
OP_CHOICES[OP_LOG] = 'log'
OP_CHOICES[OP_LOG10] = 'log10'
OP_CHOICES[OP_SQRT] = 'sqrt'
OP_CHOICES[OP_ABS] = 'abs'
OP_CHOICES[OP_POW] = 'pow'
OP_CHOICES[OP_COLUNM] = 'colunm'
OP_CHOICES[OP_SLICE] = 'slice'
OP_CHOICES[OP_ROW] = 'row'
OP_CHOICES[OP_RADIAL] = 'radial average'
OP_CHOICES[OP_RESET] = 'reset'


binaryCondition = ('(operation == %d or operation == %d or operation == %d or '
                   'operation == %d or operation == %d or operation == %d) ' %
                   (OP_PLUS, OP_MINUS, OP_MULTIPLY,
                    OP_DIVIDE, OP_MINIMUM, OP_MAXIMUM))

#noValueCondition = '(operation == 7 or operation == 8 or operation == 9 or '\
#                   'operation == 10 or operation == 15 or operation == 16) '
noValueCondition = '(operation == %d or operation == %d or operation == %d or '\
                   'operation == %d or operation == %d or operation == %d) '%\
                   (OP_LOG, OP_LOG10, OP_SQRT,\
                    OP_ABS, OP_POW, OP_RESET)

intValueCondition = '(operation == %d or operation == %d)'%(OP_COLUNM, OP_ROW)

dotCondition = 'operation == %d'%OP_DOTPRODUCT
powCondition = 'operation == %d'%OP_POW

operationDict = {OP_PLUS : ' --plus ', OP_MINUS : ' --minus ',
                 OP_MULTIPLY : ' --mult ', OP_DIVIDE : ' --divide ',
                 OP_MINIMUM : ' --min ', OP_MAXIMUM : ' --max ',
                 OP_DOTPRODUCT : ' --dot_product ', OP_LOG : ' --log ',
                 OP_LOG10 : ' --log10', OP_SQRT : ' --sqrt ',
                 OP_ABS : ' --abs ', OP_POW : ' --pow ',
                 OP_SLICE : ' --slice ',  OP_RADIAL : ' --radial_avg ',
                 OP_RESET : ' --reset ', OP_COLUNM: '--column',
                 OP_ROW: '--row'}


class XmippOperateHelper():
    """ Some image operations such as: Dot product or Summation. """
    _label = 'image operate'

    def __init__(self, **args):
        self._program = "xmipp_image_operate"

    #--------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        form.addParam('operation', params.EnumParam, choices=OP_CHOICES.keys(),
                      default=OP_PLUS,
                      label="Operation",
                      help="Binary operations: \n"
                           "*plus*: Sums two images, volumes or adds a "
                           "numerical value to an image. \n"
                           "*minus*: Subtracts two images, volumes or "
                           "subtracts a numerical value to an image. \n"
                           "*multiply*: Multiplies two images, volumes, or "
                           "multiplies per a given number. \n"
                           "*divide*: Divides two images, volumes, or divides "
                           "per a given number. \n"
                           "*minimum*: Minimum of two images, volumes, or "
                           "number (pixel-wise). \n"
                           "*maximum*: Maximum of two images, volumes, or "
                           "number (pixel-wise). \n"
                           "*dot product*: Dot product between two images or"
                           " volumes. \n"
                           "Unary operations: \n"
                           "*log*: Computes the natural logarithm of an "
                           "image. \n"
                           "*log10*: Computes the decimal logarithm of "
                           "an image. \n"
                           "*sqrt*: Computes the square root of an image \n"
                           "*abs*: Computes the absolute value of an image. \n"
                           "*pow*: Computes the power of an image. \n"
                           "*slice*: Extracts a given slice from a volume "
                           "(first slice=0). \n"
                           "*column*: Extracts a given column from a image "
                           "or volume. \n"
                           "*row*: Extracts a given row from a image or "
                           "volume. \n"
                           "*radial average*: Compute the radial average of "
                           "an image. \n"
                           "*reset*: Set the image to 0")
        form.addParam('isValue', params.BooleanParam, default=False,
                      label="Second operand is a value?",
                      condition=binaryCondition,
                      help="Set to true if you want to use a "
                           "value of the second operand")
        self._defineSpecificParams(form)
        form.addParam('value', params.FloatParam,
                      allowNull=True,
                      condition='isValue and %s or %s' %
                                (binaryCondition, powCondition),
                      label='Input value ',
                      help = 'Set the desire float value')
        form.addParam('intValue', params.IntParam,
                      allowNull=True,
                      condition=intValueCondition,
                      label='Input value ',
                      help = 'This value must be integer')
    
    #--------------------------- INSERT STEPS functions ------------------------
    def _insertProcessStep(self):
        operationStr = operationDict[self.operation.get()]
        self._insertFunctionStep("operationStep", operationStr)
    
    #--------------------------- UTILS functions ------------------------------
    def _isBinaryCond(self):
        operation = self.operation.get()
        return (operation == OP_PLUS or operation == OP_MINUS or
                operation == OP_MULTIPLY or operation == OP_DIVIDE
                or operation == OP_MINIMUM or operation == OP_MAXIMUM)
    
    def _isNoValueCond(self):
        operation = self.operation.get()
        return (operation == OP_LOG or operation == OP_LOG10 or
                operation == OP_SQRT or operation == OP_ABS or
                operation == OP_RADIAL or operation == OP_RESET)
    
    def _isPowCond(self):
        operation = self.operation.get()
        return operation == OP_POW
    
    def _isDotCond(self):
        operation = self.operation.get()
        return operation == OP_DOTPRODUCT
        
    
    def _isintValueCond(self):
        operation = self.operation.get()
        return (operation == OP_SLICE or operation == OP_COLUNM
                or operation == OP_ROW)
    
    def _getSecondSetFn(self):
        return self._getTmpPath("images_to_apply.xmd")


class XmippProtImageOperateParticles(ProtOperateParticles,
                                     XmippProcessParticles,
                                     XmippOperateHelper):
    """ Apply an operation to two sets of particles  """
    _label = 'operate particles'
    
    def __init__(self, **args):
        ProtOperateParticles.__init__(self, **args)
        XmippProcessParticles.__init__(self)
        XmippOperateHelper.__init__(self, **args)

    #--------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        XmippOperateHelper._defineProcessParams(self, form)
    
    def _defineSpecificParams(self, form):
        form.addParam('inputParticles2', params.PointerParam,
                      allowNull=True,
                      condition=(binaryCondition+' and (not isValue) or '
                                 +dotCondition),
                      label='Input Particles (2nd)',
                      help = 'Set a SetOfParticles. The particles must be '
                             'the same dimensions as the input particles.',
                      pointerClass='SetOfParticles')
    
    #--------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ convert to Xmipp image model"""
        writeSetOfParticles(self.inputParticles.get(), self.inputFn,
                            alignType=ALIGN_NONE)
        
        if self.inputParticles2.get() is not None:
            writeSetOfParticles(self.inputParticles2.get(),
                                self._getSecondSetFn(),
                                alignType=ALIGN_NONE)
    
    def operationStep(self, operationStr):
        dictImgFn = {"inputFn" : self.inputFn}
        args = self._args  % dictImgFn + operationStr
        
        if self._isBinaryCond():
            if self.isValue:
                args += ' %f' % self.value.get()
            else:
                args += ' %s' % self._getSecondSetFn()
        elif self._isPowCond():
            args += ' %f' % self.value.get()
        elif self._isDotCond():
            args += ' %s' % self._getSecondSetFn()
        elif self._isintValueCond():
            args += ' %d' % self.intValue.get()
        
        args += " -o %s --save_metadata_stack %s" % (self.outputStk,
                                                     self.outputMd)
        args += " --keep_input_columns"
        self.runJob(self._program, args)
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        def _errorDimensions():
            if not (self.inputParticles.get() is None and
                    self.inputParticles2.get() is None):
                dim1 = self.inputParticles.get().getDimensions()
                dim2 = self.inputParticles2.get().getDimensions()
                if dim1 != dim2:
                    errors.append("The number of images in two operands are "
                                  "not the same. ")
        def _errorSize():
            if not (self.inputParticles.get() is None and
                    self.inputParticles2.get() is None):
                size1 = self.inputParticles.get().getSize()
                size2 = self.inputParticles2.get().getSize()
                if size1 != size2:
                    errors.append("Size of both *SetOfParticles* are not the"
                                  " same. ")
        
        if self._isBinaryCond():
            if not self.isValue:
                if not self.inputParticles2.get() is None:
                    _errorSize()
                    _errorDimensions()
        elif self._isDotCond():
            if not self.inputParticles2.get() is None:
                _errorDimensions()
        return errors


class XmippProtImageOperateVolumes(ProtOperateVolumes,
                                   XmippProcessVolumes,
                                   XmippOperateHelper):
    """ Apply an operation to two sets of volumes """
    _label = 'operate volumes'
     
    def __init__(self, **args):
        ProtOperateVolumes.__init__(self, **args)
        XmippProcessVolumes.__init__(self)
        XmippOperateHelper.__init__(self, **args)
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        XmippOperateHelper._defineProcessParams(self, form)
    
    def _defineSpecificParams(self, form):
        form.addParam('inputVolumes2', params.PointerParam,
                      allowNull=True,
                      condition=(binaryCondition + ' and (not isValue) or '
                                 + dotCondition),
                      label='Input Volumes (2nd)',
                      help = 'This parameter depends of the input volume(s). '
                             'If it is set a volume (or a SetOfVolumes) as '
                             'input, this must be a *Volume* (or '
                             '*SetOfVolumes*) object.',
                      pointerClass='Volume, SetOfVolumes')

    #--------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ convert to Xmipp image model"""
        if not self._isSingleInput():
            writeSetOfVolumes(self.inputVolumes.get(), self.inputFn)
            
            if self.inputVolumes2.get() is not None:
                writeSetOfVolumes(self.inputVolumes2.get(),
                                  self._getSecondSetFn())
    
    def operationStep(self, operationStr):
        dictImgFn = {"inputFn" : self.inputFn}
        args = self._args  % dictImgFn + operationStr
        
        if self._isBinaryCond():
            if self.isValue:
                args += ' %f' % self.value.get()
            else:
                args += ' %s' % self._getSecondVolumeFn()
        elif self._isPowCond():
            args += ' %f' % self.value.get()
        elif self._isDotCond():
            args += ' %s' % self._getSecondVolumeFn()
        elif self._isintValueCond():
            args += ' %d' % self.intValue.get()
        
        args += " -o %s " % self.outputStk
        if not self._isSingleInput():
            args += " --save_metadata_stack %s" \
                    " --keep_input_columns" % self.outputMd
        self.runJob(self._program, args)
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        def _errorDimensions():
            if not (self.inputVolumes.get() is None and
                    self.inputVolumes2.get() is None):
                dim1 = self.inputVolumes.get().getDimensions()
                dim2 = self.inputVolumes2.get().getDimensions()
                if dim1 != dim2:
                    errors.append("The number of volumes in two operands are "
                                  "not the same. ")
        def _errorSize():
            if (not (self.inputVolumes.get() is None and
                    self.inputVolumes2.get() is None) and
                    not self._isSingleInput()):
                size1 = self.inputVolumes.get().getSize()
                size2 = self.inputVolumes2.get().getSize()
                if size1 != size2:
                    errors.append("Size of both volumes are not the same. ")
        
        if self._isBinaryCond():
            if not self.isValue:
                if not self.inputVolumes2.get() is None:
                    _errorSize()
                    _errorDimensions()
        elif self._isDotCond():
            if not self.inputVolumes2.get() is None:
                _errorDimensions()
        return errors
    
    #--------------------------- UTILS functions ------------------------------
    def _getSecondVolumeFn(self):
        if self._isSingleInput():
            return getImageLocation(self.inputVolumes2.get())
        else:
            return self._getSecondSetFn()
