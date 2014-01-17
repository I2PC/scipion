# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
"""
This module define the text used in the application.
"""

class Message():
    
    # Example Usage: 
        #       MyMessage = Message()
        #       print MyMessage.label

    # Projects Template
    LABEL_PROJECTS='Projects'
    LABEL_CREATE_PROJECT='Create Project'
    TITLE_CREATE_PROJECT='Enter the project name'
    LABEL_DELETE_PROJECT='Delete Project'
    TITLE_DELETE_PROJECT='Confirm project deletion'
    LABEL_MODIFIED='Modified: '
    
    # Project Content Template
    LABEL_PROJECT='Project '
    #-- Protocol Treeview --
    LABEL_WORKFLOW='Workflow View: '
    LABEL_PROTTREE_NONE='No athletes.'
    #-- Toolbar --
    LABEL_EDIT='Edit'
    LABEL_COPY='Copy'
    LABEL_DELETE='Delete'
    LABEL_BROWSE='Browse'
    LABEL_STOP='Stop'
    LABEL_ANALYZE='Analyze Results'
    LABEL_TREE='Tree'
    LABEL_LIST='List'
    LABEL_REFRESH='Refresh'
    #-- Tabs --
    LABEL_DATA='Data'
    LABEL_SUMMARY='Summary'
    LABEL_INPUT='Inputs'
    LABEL_OUTPUT='Outputs'
    NO_INFO_SUMMARY = 'No summary information.'
    
    # Protocol Form Template
    TITLE_NAME_RUN=' Protocol Run: '
    TITLE_RUN='Run'
    TITLE_RUN_NAME='Run name'
    LABEL_COMMENT='Describe your run here...'
    TITLE_RUN_MODE='Run mode'
    LABEL_RUN_MODE_RESUME='resume'
    LABEL_RUN_MODE_RESTART='restart'
    TITLE_EXPERT='Expert Level'
    LABEL_EXPERT_NORMAL='Normal'
    LABEL_EXPERT_ADVANCE='Advanced'
    LABEL_EXPERT_EXPERT='Expert'
    TITLE_EXEC='Execution'
    TITLE_EXEC_HOST='Execution host'
    TITLE_THREADS='Threads'
    TITLE_MPI='MPI'
    TITLE_QUEUE='Launch to queue?'
    LABEL_QUEUE_YES='Yes'
    LABEL_QUEUE_NO='No'
    LABEL_PARAM_YES='Yes'
    LABEL_PARAM_NO='No'
    LABEL_BUTTON_CLOSE='Close'
    LABEL_BUTTON_SAVE='Save'
    LABEL_BUTTON_EXEC='Execute'
    LABEL_BUTTON_VIS='Visualize'
    LABEL_BUTTON_RETURN='Save'
    
    TITLE_LAUNCHED='Success'
    LABEL_LAUNCHED='The protocol was launched successfuly'
    LABEL_FOUND_ERROR='Errors found'
    TITLE_SAVED='Success'
    LABEL_SAVED='The protocol was saved successfuly'
    
    
    
    