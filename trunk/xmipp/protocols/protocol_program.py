#!/usr/bin/env python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

#------------------------------------------------------------------------------------------------
# Generic protocol for all Xmipp programs


from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict

class ProtProgram(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.name, scriptname, project)
        self.Import = 'from protocol_ml2d import *'

    def validate(self):
        return []
        #return ["Protocol not implemented yet..."]
    
    def summary(self):
        return ["This is a test summary",
                "Need a real summary here..."
                ]
        
    def defineSteps(self):
        print '*********************************************************************'
        progId = "ml"
        if (self.DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId
#        action = "Executing"
#        if (restart):
#            action = "Restarting"
#            
#        print '*  %s %s program :' % (action, program)
        restart = False
        if (restart):
            pass 
            #Not yet implemented
            #params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            params = ' -i %s --oroot %s/%s2d' % (self.InSelFile, self.WorkingDir, progId)
            # Number of references will be ignored if -ref is passed as expert option
            if self.ExtraParams.find("--ref") == -1:
                params += ' --nref %i' % self.NumberOfReferences
            params += ' ' + self.ExtraParams
            if (self.DoFast and not self.DoMlf):
                params += ' --fast'
            if (self.DoNorm):
                params += ' --norm'
            if (self.DoMirror):
                params += ' --mirror'
            if (self.NumberOfThreads > 1  and not self.DoMlf):
                params += ' --thr %i' % self.NumberOfThreads
            if (self.DoMlf):
                if (self.DoCorrectAmplitudes):
                    params += ' --ctfdat %s' % self.InCtfDatFile
                else:
                    params += ' --no_ctf --pixel_size %f' % self.PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' % self.HighResLimit
                    
        self.Db.insertStep('runJob', 
                             programname=program, 
                             params=params,
                             NumberOfMpi = self.NumberOfMpi,
                             NumberOfThreads = self.NumberOfThreads)

