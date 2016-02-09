# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), Sep 2013
# * Ported to Scipion:
# *           Vahid Abrishami (vabrishami@cnb.csic.es), Oct 2014
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

import os
import xmipp
from pyworkflow.protocol.params import FileParam, EnumParam, PointerParam
from pyworkflow.em import *
from protocol_process import XmippProcessParticles, XmippProcessVolumes
from pyworkflow.utils.properties import Message
from ..convert import writeSetOfParticles, writeSetOfVolumes

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
OP_COLUMN = 13
OP_ROW = 14
OP_RADIAL = 15
OP_AVERAGE = 16
OP_RESET = 17

OP_CHOICES = ['plus', 'minus', 'multiply', 'divide', 'minimum', 'maximum',
              'dot product', 'log', 'log10', 'sqrt', 'abs', 'pow', 'slice',
              'column', 'row', 'radial', 'average', 'reset']

conditionStr = '(operation == 0 or operation == 1 or operation == 2 or ' \
                'operation == 3 or operation == 4 or operation == 5 or ' \
                'operation == 6)'

operationDict = {OP_PLUS : ' --plus ', OP_MINUS : ' --minus ',
                 OP_MULTIPLY : ' --mult ', OP_DIVIDE : ' --divide ',
                 OP_MINIMUM : ' --min ', OP_MAXIMUM : ' --max ',
                 OP_DOTPRODUCT : ' --dot_product ', OP_LOG : ' --log ',
                 OP_LOG10 : ' --log10', OP_SQRT : ' --sqrt ',
                 OP_ABS : ' --abs ', OP_POW : ' --pow ',
                 OP_SLICE : ' --slice ', OP_COLUMN : ' --column ',
                 OP_ROW : ' --row ', OP_RADIAL : ' --radial_avg ',
                 OP_RESET : ' --reset '}


class XmippProtImageOperate():
    """ Some image operations such as: Dot product or Summation. """
    _label = 'image operate'

    def __init__(self,_isParticle, **args):
        self._program = "xmipp_image_operate"
        self._isParticle = _isParticle

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        form.addParam('operation', EnumParam, choices=OP_CHOICES,
                      default=OP_PLUS,
                      label="Operation")
        form.addParam('isValue', BooleanParam, default=False,
                      label="Second operand is a value", help="Set to true if you want to use a "
                                                              "value of the second operand")
        if self._isParticle:
            form.addParam('inputParticles2', PointerParam, important=True,
                           condition=conditionStr + ' and (not isValue)',
                           label='Input Particles (2nd)',
                           help = 'It can be a metadata, a stack of images',
                           pointerClass='SetOfParticles')
        else:
            form.addParam('inputVolumes2', PointerParam, important=True,
                           condition=conditionStr + ' and (not isValue)',
                           label='Input Volumes (2nd)',
                           help = 'It can be a metadata, a stack of volumes',
                           pointerClass='Volume, SetOfVolumes')
        form.addParam('inputValue', FloatParam, important=True,
                       condition='isValue', label='Input value ',
                       help = 'For the operations which support values')


    def _insertProcessStep(self):

        inputFn = self.inputFn
        # determine the command line for
        operationStr = operationDict[self.operation.get()]
        if not self.isValue.get():
            if self._isParticle:
                inputFn2 = self._getTmpPath('input_particles2.xmd')
                writeSetOfParticles(self.inputParticles2.get(), inputFn2)
            else:
                inputFn2 = self._getTmpPath('input_volumes2.xmd')
                #if we do not create this case them getAlignment must be defined
                if isinstance(self.inputVolumes2.get(), Volume):
                    inputFn2 = self.inputVolumes2.get().getFileName()
                else:#set of volumes
                    #TODO ROB: how to deal with binary conversions?
                    writeSetOfVolumes(self.inputVolumes2.get(), inputFn2)
            args = (self._args + operationStr + ' ' + inputFn2) % locals()
        else:
            val = self.inputValue.get()
            args = (self._args + operationStr) % locals()
            args += ' %f' % val
        self._insertFunctionStep("operationStep", args)
        #--------------------------- INSERT steps functions --------------------------------------------


class XmippProtImageOperateParticles(ProtOperateParticles, XmippProcessParticles, XmippProtImageOperate):
    """ Apply an operation to two sets of particles  """
    _label = 'calculator2D'
    _isParticle = True

    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
        XmippProcessParticles.__init__(self)
        XmippProtImageOperate.__init__(self, self._isParticle, **args)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtImageOperate._defineProcessParams(self, form)

    #--------------------------- STEPS functions ---------------------------------------------------
    def operationStep(self, args):
        print "operationStep"
        args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)
        self.runJob("xmipp_image_operate", args)
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        operation = self.operation.get()
        N1 = self.inputParticles.get().getSize()
        if self.isValue.get == False:
            N2 = self.inputParticles2.get().getSize()
        else:
            N2=1
        checkDimension = False
        if operation == OP_COLUMN or operation == OP_SLICE or operation == OP_ROW:
            if not self.isValue.get():
                errors.append("You should give a number for the column, slice or row")
        elif operation == OP_DOTPRODUCT:
            if self.isValue.get():
                errors.append("Second operand can not be a number")
            else:
                checkDimension = True
        elif operation == OP_PLUS or operation == OP_MINUS or operation == OP_MULTIPLY or \
                          operation == OP_DIVIDE or operation == OP_MINIMUM or \
                          operation == OP_MAXIMUM:
            if not self.isValue.get():
                checkDimension = True
        if checkDimension:
            x1, y1, z1 = self.inputParticles.get().getDimensions()
            x2, y2, z2 = self.inputParticles2.get().getDimensions()
            if x1 != x2 or y1 != y2 or z1 != z2:
                errors.append("Image sizes in the two operands are not the same")
            if N2 > 1:
                if N1 != N2:
                    errors.append("The number of images in two operands are not the same")
        return errors


class XmippProtImageOperateVolumes(ProtOperateVolumes,
                                   XmippProcessVolumes,
                                   XmippProtImageOperate):
    """ Apply an operation to two sets of volumes """
    _label = 'calculator3D'
    _isParticle = False

    def __init__(self, **args):
        ProtPreprocessVolumes.__init__(self, **args)
        XmippProcessVolumes.__init__(self)
        XmippProtImageOperate.__init__(self, self._isParticle, **args)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        XmippProtImageOperate._defineProcessParams(self, form)

    #--------------------------- STEPS functions ---------------------------------------------------
    def operationStep(self, args):
       if self._isSingleInput():
            args += " -o %s" % ImageHandler().getVolFileName(self.outputStk)
       else:
            args += " -o %s --save_metadata_stack %s --keep_input_columns" % (self.outputStk, self.outputMd)

       self.runJob("xmipp_image_operate", args)
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        operation = self.operation.get()
        #N1 = self.inputVolumes.get().getSize()
        #N2 = self.inputVolumes2.get().getSize()
        checkDimension = False
        if operation == OP_COLUMN or operation == OP_SLICE or operation == OP_ROW:
            if not self.isValue.get():
                errors.append("You should give a number for the column, slice or row")
        elif operation == OP_DOTPRODUCT:
            if self.isValue.get():
                errors.append("Second operand can not be a number")
            else:
                checkDimension = True
        elif operation == OP_PLUS or operation == OP_MINUS or operation == OP_MULTIPLY or \
                          operation == OP_DIVIDE or operation == OP_MINIMUM or \
                          operation == OP_MAXIMUM:
            if not self.isValue.get():
                checkDimension = True
        if checkDimension:
            #TODO ROB: this should be getDimensions even for single objects
            #Do the same with 2D?
            # but I cannot make it work
            if isinstance(self.inputVolumes2.get(), Volume)\
               or isinstance(self.inputVolumes2.get(), Particle):
                x1, y1, z1 = self.inputVolumes.get().getDim()
                x2, y2, z2 = self.inputVolumes2.get().getDim()
            else:
                x1, y1, z1 = self.inputVolumes.get().getDimensions()
                x2, y2, z2 = self.inputVolumes2.get().getDimensions()
            if x1 != x2 or y1 != y2 or z1 != z2:
                errors.append("Volume sizes in the two operands are not the same")
            #if N2 > 1:
                #if N1 != N2:
                    #errors.append("The number of volumes in two operands are not the same")
        return errors
