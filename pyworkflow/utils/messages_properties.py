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
    label_projects='Projects'
    label_create_project='Create Project'
    title_create_project='Enter the project name'
    label_delete_project='Delete Project'
    title_delete_project='Confirm project deletion'
    label_modified='Modified: '
    
    # Project Content Template
    label_project='Project '
    #-- Protocol Treeview --
    label_workflow='Workflow View: '
    label_protTree_None='No athletes.'
    #-- Toolbar --
    label_edit='Edit'
    label_copy='Copy'
    label_delete='Delete'
    label_browse='Browse'
    label_stop='Stop'
    label_analyze='Analyze Results'
    label_tree='Tree'
    label_list='List'
    label_refresh='Refresh'
    #-- Tabs --
    label_data='Data'
    label_summary='Summary'
    label_input='Inputs'
    label_output='Outputs'
    no_info_summary = 'No summary information.'
    
    # Protocol Form Template
    title_name_run=' Protocol Run: '
    title_run='Run'
    title_run_name='Run name'
    label_comment='Describe your run here...'
    title_run_mode='Run mode'
    label_run_mode_resume='resume'
    label_run_mode_restart='restart'
    title_expert='Expert Level'
    label_expert_normal='Normal'
    label_expert_advance='Advanced'
    label_expert_expert='Expert'
    title_exec='Execution'
    title_exec_host='Execution host'
    title_threads='Threads'
    title_mpi='MPI'
    title_queue='Launch to queue?'
    label_queue_yes='Yes'
    label_queue_no='No'
    label_param_yes='Yes'
    label_param_no='No'
    label_button_close='Close'
    label_button_save='Save'
    label_button_exec='Execute'
    label_button_vis='Visualize'
    label_button_return='Save'
    
    title_launched='Success'
    label_launched='The protocol was launched successfuly'
    label_found_error='Errors found'
    title_saved='Success'
    label_saved='The protocol was saved successfuly'
    
    
    
    