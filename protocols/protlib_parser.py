#!/usr/bin/env xmipp_python
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
 ***************************************************************************

This module contains the classes to perform the protocol header
parsing, getting protocol variables and meta-information. 
 '''

import os, sys
from protlib_include import *
from protlib_xmipp import greenStr, redStr

# Variable types
VAR_STRING = 0
VAR_TEXT = 1 # Multiline string
VAR_BOOL = 2
VAR_NUMBER = 3 # Int or float
        
class ProtocolVariable():
    '''Store meta-information about variables.'''
    def __init__(self, parser):
        self.name = None
        self.value = None
        self.comment = None
        self.commentline = None
        self.help = None
        self.tags = {}
        self.conditions = {} #map with conditions to be visible
        self.validators = [] #list with validators for this variable
        self.type = None
        self.section = None # Section to which belongs the variable
         
        self.parser = parser
        self.childs = [] # Only section will use childs    
                  
        self.tkvar = None #used from gui
        self.tktext = None #special text widget     
    
    def setTags(self, tags):
        for k, v in tags:
            self.tags[k] = v
            
    def hasTag(self, tag):
        return tag in self.tags
    
    def setTag(self, tagKey, tagValue=None):
        self.tags[tagKey] = tagValue
        
    def getTag(self, tagKey):
        return self.tags[tagKey]
    
    def getTagValues(self, tagKey):
        return [v.strip() for v in self.tags[tagKey].split(',')]
    
    def getValue(self):
        #print "%s " % redStr(self.name)
        #print " = %s" % greenStr(self.value)
        return self.value
    
    def getTkValue(self):
        '''this function will be used with gui,
        it will update variable value with the
        gui tkvar value '''
        if self.isHidden():
            return self.value
        if self.tkvar:
            return self.tkvar.get()
        elif self.tktext:
            return self.tktext.get(1.0, 'end')
        else:
            return self.value
        
    def updateValue(self):
        self.value = self.getTkValue()
            
    def setTkValue(self, value):
        '''Update the this function will be used with gui,
        it will update variable value with the
        gui tkvar value '''
        if self.tkvar:
            self.tkvar.set(value)
        elif self.tktext:
            self.tktext.clear()
            self.tktext.addText(value)
        else:
            self.value = value
    
    def setValue(self, value):
        self.value = value
    
    def getLiteralValue(self):
        if self.isText():
            return '"""%s"""' % self.getValue()
        elif self.isString():
            return '"%s"' % self.getValue()
        return self.getValue()

    def checkType(self):
        '''Check the type of the variable depending on
        the literal value read from header '''
        v = startswithQuotes(self.value, 3)
        if not v is None:
            self.type = VAR_TEXT
            self.value = v
            return
        v = startswithQuotes(self.value, 1)
        if not v is None:
            self.type = VAR_STRING
            self.value = v
            return
        if self.value in ['True', 'False']:
            self.type = VAR_BOOL
            return 
        if isNumber(self.value):
            self.type = VAR_NUMBER
            return
        raise Exception('Invalid value: "%s" for variable: "%s"' % (self.value, self.name))
    
    def isExpert(self):
        return self.hasTag('expert')
    
    def isVisualize(self):
        if self.section is None:
            return self.hasTag('visualize')
        else:
            return self.section.hasTag('visualize')
    
    def isSection(self):
        return self.hasTag('section')
    
    def isHidden(self):
        return self.hasTag('hidden') or self.isCite()
    
    def isCite(self):
        return self.hasTag('cite')
    
    def isRun(self):
        return self.hasTag('run')
    
    def isString(self):
        return self.type == VAR_STRING
    
    def isText(self):
        return self.type == VAR_TEXT
    
    def isBoolean(self):
        return self.type == VAR_BOOL
    
    def isNumber(self):
        return self.type == VAR_NUMBER
    
    def hasQuestion(self):
        return self.hasTag('has_question')
    
    def hasValidate(self):
        return self.hasTag('validate')
    
    def isSeparator(self):
        return self.commentline.startswith("#---")

    def getLines(self):   
        self.updateValue()
        lines = ['']
        if self.isSection():
            sepLine = '#%s' % ('-' * 100) 
            lines.append("%s\n%s\n%s" % (sepLine, self.commentline, sepLine))
            for v in self.childs:
                lines += v.getLines()
        else:    
            lines.append(self.commentline)
            if self.help:
                lines.append(self.help)
            if self.isText(): #Store the value in the comment field
                template = '%s = """\n%s\n"""'            
            elif self.isString():
                template = '%s = "%s"\n'
            else:
                template = '%s = %s\n'
            lines.append(template % (self.name, self.getValue()))
        return lines
    
    def addChild(self, variable):
        '''Add child variable to section '''
        if self.isSection():
            self.childs.append(variable)
            variable.section = self
        else:
            raise Exception("Adding child variable to non-section variable")
    
    def getSection(self):
        return self.section
    
    def satisfiesCondition(self):
        self.updateValue()
        if not self.hasTag('condition'):
            return True
        import re
        condition = self.tags['condition']
        tokens = re.split('\W+', condition)
        for t in tokens:
            if self.parser.hasVariable(t):
                condition = condition.replace(t, self.parser.getVariable(t).getLiteralValue())
                condition = condition.replace('[', '(')
                condition = condition.replace(']', ')')
        if len(condition):#,  eval(condition)
            return eval(condition)
        return True 
    
    def hasValidator(self, validator):
        '''check if a validator is already present'''
        return 'validator' + validator in self.validators
    
    def addValidator(self, validator):
        ''' Add a validator to the list if not present.
        The name should be without the 'validator' prefix, that will 
        be added automatically '''
        if not self.hasValidator(validator):
            self.validators.append('validator' + validator) 
        
    def updateValidators(self):
        ''' Take validators names comming in 'validate' tag and 
        add the 'validator' prefix, storing all validators names
        in a list'''
        if self.hasValidate():
            for v in self.getTagValues('validate'):
                self.addValidator(v)
        if self.isNumber() and not self.hasValidator('IsInt'):
            self.addValidator('IsFloat')
        if self.isRun():
            self.addValidator('ValidRun')
            
            
class ProtocolParser():
    ''' Class to parse the protocol header files and extract the
    variables information arranged in sections '''
    def __init__(self, header_script):
        self.script = header_script
        self.readProtocolScript(header_script)
        self.parseHeaderLines()

    def readProtocolScript(self, script):
        ''' Read header lines and split in three main categories:
        before header
        header
        after header'''
        self.pre_header_lines = []
        self.header_lines = []
        self.post_header_lines = []
        begin_of_header = False
        end_of_header = False   
        has_expert = False
        
        try:   
            f = open(script, 'r')
        except Exception, e:
            raise Exception("Script read failed", "Couldn't read from script file <%s>" % script)
        
        for line in f:
            if '{begin_of_header}' in line:
                begin_of_header = True
            elif '{end_of_header}' in line:
                if not has_expert:
                    self.header_lines += expandExpert().splitlines()
                end_of_header = True
            else:
                if not begin_of_header:
                    self.pre_header_lines.append(line)
                elif not end_of_header:
                    if 'ShowExpertOptions' in line:
                        has_expert = True                    
                    if '{eval}' in line:
                        evalStr = line.split('{eval}')[1].strip()
                        linesStr = eval(evalStr)
                        if linesStr:
                            self.header_lines += linesStr.splitlines()                    
                    else:
                        self.header_lines.append(line)
                else:
                    self.post_header_lines.append(line)                
        f.close()
        
        if not begin_of_header:
            raise Exception('{begin_of_header} tag not found in protocol script: %s' % script)
        if not end_of_header:
            raise Exception('{end_of_header} tag not found in protocol script: %s' % script)         
            
    def moreLines(self):
        '''Check the index has not reach the end '''
        return self.index < self.count
        
    def currentLine(self):
        ''' Return the line currently parsed, using the index'''
        return self.header_lines[self.index].strip()
    
    def getLine(self):
        '''Get current line and move'''
        line = self.currentLine()
        self.nextLine()
        return line
    
    def nextLine(self, deltha=1):
        '''Move the index to next line '''
        self.index += deltha
        
    def parseMultilineString(self, isComment=True):
        '''Parse a multiline python string literal'''
        strValue = ''
        if self.moreLines():
            if isComment:
                parse = self.currentLine().startswith('"""')
            else:
                parse = '"""' in self.currentLine()
            if parse:
                pos = self.currentLine().find('"""')
                line = self.getLine()[pos:]
                condition = len(line) > 4
                while self.moreLines():
                    if len(line) > 0:
                        strValue += line + '\n'
                    if condition and line.endswith('"""'):
                        break
                    line = self.getLine()
                    condition = True
        return strValue
           
    def addSection(self, v):
        self.sections.append(v)
        self.lastSection = v
        
    def addVariable(self, v):
        self.variables[v.name] = v
        self.lastSection.addChild(v)
        
    def hasVariable(self, varName):
        return varName in self.variables
    
    def getVariable(self, varName):
        return self.variables[varName]
    
    def getValue(self, varName):
        return self.variables[varName].getValue()
    
    def getTkValue(self, varName):
        return self.variables[varName].getTkValue()
    
    def getTkBoolValue(self, varName):
        return self.getTkValue(varName) == 'True'
    
    def getLiteralValue(self, varName):
        return self.variables[varName].getLiteralValue()
    
    def setValue(self, varName, value):
        self.variables[varName].setValue(value)
        
    def parseHeaderLines(self):
        '''Parse header lines storing the information
        in the variables dictionary '''
        self.variables = {}
        self.sections = [] 
        self.lastSection = None

        import re
        #Comment regex, match lines starting by # and followed by tags with values
        #of the form {tag}(value) and ending with the comment for the GUI    
        reComment = re.compile('#\s*((?:{\s*\w+\s*}\s*(?:\([^)]*\)\s*)?)*)?\s*(.*)')
        #This take all tags and values from previous one
        reTags = re.compile('(?:{\s*(\w+)\s*}\s*(?:\(([^)]*)\))?\s*)')
        #This is the regular expression of a Variable
        #possible values are: True, False, String with single and double quotes and a number(int or float) 
        reVariable = re.compile('(\w+)\s*=\s*("""|True|False|".*"|\'.*\'|[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?|)?')
        
        self.widgetslist = []
        self.index = 0
        self.count = len(self.header_lines)
                 
        while self.moreLines():
            line = self.getLine()
            if not line.startswith("#---"): # skip separator lines
                match = reComment.search(line) #Parse the comment line
                if match:
                    v = ProtocolVariable(self)
                    v.comment = match.group(2)
                    v.commentline = line
                    
                    if match.group(1) != '':
                        tags = reTags.findall(match.group(1))
                        v.setTags(tags)
    

                    if not v.isSection():
                        v.help = self.parseMultilineString()
                        
                        if self.moreLines():
                            line = self.currentLine()
                            match2 = reVariable.match(line)
                            if match2:
                                v.name, v.value = (match2.group(1).strip(), match2.group(2).strip())
                                #print greenStr("v.name = %s" % v.name)
                                #print redStr("v.value = %s" % v.value)
                                self.addVariable(v)
                                #self.nextLine()
                                if v.value.startswith('"""'): # Parse multiline variable
                                    v.value = self.parseMultilineString(False)
                                    if emptyString(v.value):
                                        raise Exception("Parse variable failed", 
                                                        "Expecting multiline value for variable '%s'" % v.name)
                                v.checkType()
                    else: self.addSection(v)
                        
    def writeLines(self, f, linesList):
        for line in linesList:
            print >> f, line.rstrip()
            
    def write(self, f=sys.stdout):  
        self.writeLines(f, self.pre_header_lines + ['# {begin_of_header}', ''])      
        for s in self.sections:
            self.writeLines(f, s.getLines())
        self.writeLines(f, ['', '# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE'] + self.post_header_lines)
        
    def save(self, filename=None):
        if filename is not None:
            self.script = filename
        f = open(self.script, 'w+')
        self.write(f)
        f.close()
        os.chmod(self.script, 0755)
        


#---------------- Some utilities functions --------------------------

def emptyString(value):
    return len(value.strip()) == 0

def startswithQuotes(s, n):
    '''Check if a string starting with quotes or doublequotes'''
    p = '"' * n
    if s.startswith(p):
        return s.replace(p, '')
    p = "'" * n
    if s.startswith(p):
        return s.replace(p, '')
    
    return None

def isNumber(s):
    try:
        float(s)
        return True
    except Exception, e:
        return False


if __name__ == '__main__':
    script = sys.argv[1]
    parser = ProtocolParser(script)
    #parser.write()
#    parser.setValue("Comment", """5sdfsdfff...44444jjjj 
#    This is a t
#    very very long....
#    test")'""")
#    print "Behavior1", parser.getValue("Behavior")
#    parser.setValue('Behavior', 'Restart')
#    print "Behavior2", parser.getValue("Behavior")
#    parser.setValue('DirMicrographs', 'InputKK')
    parser.save("kk.py")

    print "Finished"