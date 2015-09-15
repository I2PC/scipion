import sys

def returnOwner( value):
    
    if (value == "rmarabini"):
        return "RM"
    elif (value == "coss"):
        return "COSS"
    else:
        return "RM"

def parse_File( inFileName, outFileName):
    
    # Test cases dictionary.
    test_Cases = {}
    mpi_Parsed = {}
    
    # Open files.
    inf = open(inFileName, 'r')
    outf = open(outFileName, 'w')
        
    # Skip header.
    finished = False
    while (not finished):

        newLine = inf.readline()
        if (newLine.strip() == "<XMIPP_TESTS>"):
            finished = True
        #sys.stdin.read(1)
        
    endOfFile = False
    while (not endOfFile):

        # Get new line.
        line = inf.readline().strip()

        # Check end of file.        
        if (line == "</XMIPP_TESTS>"):
            endOfFile = True
        # Check if current line is a new program.
        elif (not (line.find("<PROGRAM name=") == -1)):
            
            fullLine = False 
            newLine = line
            
            # Read full program line.
            while (not fullLine):
          
                if (newLine.find(">") == -1):
                    newLine = inf.readline().strip()
                    line = line + " " + newLine
                else:
                    fullLine = True
            
            # Get name, mpi and owner.
            firstIndex = line.index('name="') + 6
            subNewLine = line[firstIndex:len(line)]
            lastIndex = subNewLine.index('"') + firstIndex
            name = line[firstIndex:lastIndex].strip()
            
            # Check if it s aMPI programg.
            line = line[lastIndex:len(line)].strip()
            if (line.find('mpi="') != -1):
                if (line.find('mpi="TRUE') != -1):
                    is_MPI_Program = True
                else:
                    is_MPI_Program = False
                firstIndex = line.index('mpi="') + 5
                subNewLine = line[firstIndex:len(line)]
                lastIndex = subNewLine.index('"') + firstIndex
            else:
                is_MPI_Program = False
            
            # Get program owner.
            line = line[lastIndex:len(line)].strip()
            if (line.find('owner="') != -1):
                firstIndex = line.index('owner="') + 7
                subNewLine = line[firstIndex:len(line)]
                lastIndex = subNewLine.index('"') + firstIndex
                owner = line[firstIndex:lastIndex].strip()
            
            # Get class name.
            if is_MPI_Program:
                firstIndex = name.index('xmipp_mpi_') + 10
                classNameRaw = name[firstIndex:len(name)]
            
            else:
                firstIndex = name.index('xmipp_') + 6
                classNameRaw = name[firstIndex:len(name)]
            
            # Convert input class name to formatted class name.
            fullLine = False
            className = ""
            while (not fullLine):
                
                if (classNameRaw.find("_") == -1):
                    fullLine = True
                    className = className + classNameRaw[0:1].upper() + classNameRaw[1:len(classNameRaw)]
                else:
                    firstIndex = classNameRaw.index('_')
                    className = className + classNameRaw[0:1].upper() + classNameRaw[1:firstIndex]
                    classNameRaw = classNameRaw[(firstIndex+1):len(classNameRaw)]
                    
            # Add "Mpi" string to class name if it is an MPI program.
            if is_MPI_Program:
                outf.write("class " + className + "Mpi(" + className + "):\n")
                mpi_Parsed[className] = True
            else:
                if (className in mpi_Parsed.keys()):
                    print("Class name " + className + " declared after a MPI program based on it\n")
                outf.write("class " + className + "(XmippProgramTest):\n")    
            
            # Write class header.
            outf.write("\t_owner = " + returnOwner(owner) + "\n")
            outf.write("\t@classmethod\n")
            outf.write("\tdef getProgram(cls):\n")
            outf.write("\t\treturn '" + name + "'\n\n")
            
            # Parse test cases.
            finished = False
            if not is_MPI_Program:
                testCounter = 1
            else:
                if className in test_Cases.keys(): 
                    testCounter = test_Cases[className]
                else:
                    testCounter = 1
                    
                    finished = False
                    while (not finished):
                        newLine = inf.readline()
                        if (newLine.strip() == "</PROGRAM>"):
                            finished = True
                
            while (not finished):
        
                # Read until end of program.
                newLine = inf.readline()
                if (newLine.strip() == "</PROGRAM>"):
                    finished = True
                    
                    if not is_MPI_Program:
                        test_Cases[className] = testCounter
                        
                elif (not (newLine.find("<CASE") == -1)):

                    # Read full CASE line.
                    line = newLine
                    fullLine = False
                    while (not fullLine):
                  
                        if (newLine.find("</CASE>") == -1):
                            newLine = inf.readline().strip()
                            line = line + " " + newLine
                        else:
                            fullLine = True
                            
                    line = line.replace( "\n", "")

                    # Print test case header.
                    outf.write("\tdef test_case" + str(testCounter) + "(self):\n")
                    testCounter = testCounter + 1
                    
                    # Get arguments.
                    firstIndex = line.index('arguments=') + 11
                    subNewLine = line[firstIndex:len(line)]
                    lastIndex = subNewLine.index('"') + firstIndex
                    arguments = line[firstIndex:lastIndex]
                    line = line[(lastIndex+1):len(line)]
                    
                    # Print test case arguments.
                    outf.write("\t\tself.runCase(\"" + arguments + "\",\n")
                    
                    # Check if test case requires PRERUN
                    if (not line.find("<PRERUN ") == -1):
                        command = ""
                        outf.write("\t\t\t\tprerun=[")
                        
                        firstIndex = line.index('command="') + 9
                        subNewLine = line[firstIndex:len(line)]
                        lastIndex = subNewLine.index('"') + firstIndex
                        command = line[firstIndex:lastIndex]
                        outf.write("\"" + command + "\"")
                        line = line[(lastIndex+1):len(line)]
                        outf.write("],\n")
                    
                    # Get output files.
                    parsed = False
                    firstOutput = True
                    outf.write("\t\t\t\toutputs=[")
                    while (not parsed):
                        
                        if (line.find("<FILE ") == -1):
                            parsed = True
                        else:
                                            
                            firstIndex = line.index('<FILE filename="') + 16
                            subNewLine = line[firstIndex:len(line)]
                            lastIndex = subNewLine.index('"') + firstIndex
                            output = line[firstIndex:lastIndex]
                            if (firstOutput):
                                outf.write("\"" + output + "\"")
                                firstOutput = False
                            else:
                                outf.write(",\"" + output + "\"")
                            line = line[(lastIndex+1):len(line)]
                            
                    outf.write("])\n")
            
            outf.write("\n\n")
            #sys.stdin.read(1)
                  
    # Close files.
    inf.close()
    outf.close()

#def parse_Program():
parse_File("test.xml", "test_programs_xmipp.py")