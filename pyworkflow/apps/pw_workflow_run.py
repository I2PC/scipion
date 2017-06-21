#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This module is responsible for running workflows headlessly from a json description.
"""

from pyworkflow.manager import Manager
import time

from optparse import OptionParser
from pyworkflow.object import Boolean

import pyworkflow.utils as pwutils
import json

def run(project_name, workflow, launch_timeout, location):
    manager = Manager()
    pwutils.path.makePath(manager.PROJECTS)
    project = manager.createProject(project_name, location=location)
    protocols = project.loadProtocols(workflow)

    labels_for_queue = find_protocols_to_queue(workflow)
    streaming = workflow_is_streaming(workflow)

    graph = project.getGraphFromRuns(protocols.values())
    nodes = graph.getRoot().iterChildsBreadth()

    protocols_by_id = {protocol.strId():protocol for protocol in protocols.values()}

    for node in nodes:
        name = node.getName()
        protocol = protocols_by_id[name]
        parents = collect_parents(node, protocols_by_id)

        # the id's don't get reused in the new project so comparing labels
        if protocol.getObjLabel() in labels_for_queue:
            protocol._useQueue = Boolean(True)

        launch_when_ready(parents, project, protocol, name, launch_timeout, streaming)


def find_protocols_to_queue(workflow):
    with open(workflow) as f:
        protocols = json.load(f)
        labels_for_queue = {protocol['object.label'] for protocol in protocols if
                     'useQueue' in protocol and protocol['useQueue']}
        return labels_for_queue


def workflow_is_streaming(workflow):
    with open(workflow) as f:
        protocols = json.load(f)
        return any(['dataStreaming' in protocol and protocol['dataStreaming'] for protocol in protocols])


def collect_parents(node, protocols_by_id):
    return {protocols_by_id[parent_node.getName()]
            for parent_node in node.getParents() if not parent_node.isRoot()}


def launch_when_ready(parents, project, protocol, name, launch_timeout, streaming):
    start = time.time()
    while time.time() < start + launch_timeout:
        update_protocols(parents, project)

        errors = protocol.validate()
        parents_not_finished = any([not parent.isFinished() for parent in parents])

        if errors:
            print "waiting to launch protocol " + name
            print "current validation errors are: "
            print errors
        elif not streaming and parents_not_finished:
            print "waiting to launch protocol " + name + " as not all parents finished"
            print "statuses are " + str([parent.getStatus() for parent in parents])
        else:
            print "launching protocol " + name
            project.launchProtocol(protocol)
            break

        check_protocol(protocol)
        for parent in parents:
            check_protocol(parent)
        time.sleep(1)

def update_protocols(parents, project):
    for parent in parents:
        project._updateProtocol(parent)


def check_protocol(protocol):
    if protocol.isFailed():
        print "\n>>> ERROR running protocol %s" % protocol.getRunName()
        print "    FAILED with error: %s\n" % protocol.getErrorMessage()
        raise Exception("ERROR launching protocol.")


if __name__ == '__main__':
    parser = OptionParser(usage='scipion workflow -p <project name> -w <workflow> -l <location of project>')
    parser.add_option("-p", "--project", dest="project_name", help="the name of the new project")
    parser.add_option("-w", "--workflow", dest="workflow", help="the path of the workflow json")
    parser.add_option("-t", "--timeout", dest="timeout", help="timeout for launching protocols", default=60 * 60)
    parser.add_option("-l", "--location", dest="location", help="location of project")

    (opt, args) = parser.parse_args()

    if None in [opt.project_name, opt.workflow, opt.timeout, opt.location]:
        parser.error('need to specify name, workflow, and location')

    run(opt.project_name, opt.workflow, opt.timeout, opt.location)
