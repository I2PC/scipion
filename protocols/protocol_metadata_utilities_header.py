#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for performing operations over metadatas 
#
# Author: Javier Vargas, Sept 2013 
#
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------

# {file}(*.xmd){validate}(PathExists) Input metadata to process:
""" 
Provide a metadata file
"""
MetadataInput = ""

""" 
Different types of operations on metadata files
"""
#  {list_combo} (set, operate, query, fill) Select an option
OperationType = "set"
# {list_combo}(union,union_all,intersection,subtraction,join,natural_join,merge){condition}(OperationType=="set") Set Operation
OperationTypeSet = "union"

""" 
Metadata to operate on the input metada file
"""
# {file}{condition}(OperationType =="set") Metadata file to operate
md2 = ""

# {condition}(OperationType =="set") Label
Label = ""

# {list_combo}(sort,random_subset,bootstrap,randomize,keep_column,drop_column,rename_column,modify_values){condition}(OperationType=="operate") Operation
Operation_TypeOp = "sort"

# {condition}(OperationType =="operate" and Operation_TypeOp=="sort") Sort label
LabelOp = ""

# {condition}(OperationType =="operate" and Operation_TypeOp=="random_subset") Size
LabelSize = ""

# {condition}(OperationType =="operate" and Operation_TypeOp=="keep_column") Keep some columns using label
Label_KC = ""

# {condition}(OperationType =="operate" and Operation_TypeOp=="drop_column") Drop some columns using label
Label_DC = ""

# {condition}(OperationType =="operate" and Operation_TypeOp=="rename_column") Rename a column
Label_RC = ""

# {condition}(OperationType =="operate" and Operation_TypeOp=="modify_values") SQLite expression to modify the metadata
Expression = ""

# {list_combo}(copy,move,delete,import_txt){condition}(OperationType=="file")  File operations
File_Operation = "copy"

# {condition}(OperationType =="file" and File_Operation=="copy") Copy files in metadata md1 to directory path
DirectoryCopy = ""
# {condition}(OperationType =="file" and File_Operation=="copy") Copy files at label column
LabelCopy = ""

# {condition}(OperationType =="file" and File_Operation=="move") Move files in metadata md1 to directory path
DirectoryMove = ""
# {condition}(OperationType =="file" and File_Operation=="move") Move files at label column
LabelMove = ""

# {condition}(OperationType =="file" and File_Operation=="delete") Delete files at label column
LabelDirectory = ""

# {condition}(OperationType =="file" and File_Operation=="import_txt") Import a text file specifying its columns
LabelImport = ""

# {list_combo}(select,count,sum){condition}(OperationType=="query")  Query operations
Query_Operation = "select"
# {condition}(OperationType =="query" and Query_Operation=="select") New metadata satisfying expression
SelectExpression = ""
# {condition}(OperationType =="query" and Query_Operation=="count") Label
LabelCount = ""
# {condition}(OperationType =="query" and Query_Operation=="sum") Label 1
LabelSum1 = ""
# {condition}(OperationType =="query" and Query_Operation=="sum") Label 2
LabelSum2 = ""

# {list_combo}(constant,lineal,rand_uniform,rand_gaussian,rand_student){condition}(OperationType=="fill")  Fill operations
Fill_Operation = "constant"
# {condition}(OperationType =="fill") Label
LabelFill = ""
# {condition}(OperationType =="fill" and Fill_Operation=="constant") Value
ValueCont = ""
# {condition}(OperationType =="fill" and Fill_Operation=="lineal") Init value
InitValueLinealCont1 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="lineal") Step
InitValueLinealCont2 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_uniform") A value
InitValueRandUnifCont1 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_uniform") B value
InitValueRandUnifCont2 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_gaussian") Mean value
InitValueRandGaussCont1 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_gaussian") Stddev value
InitValueRandGaussCont2 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_student") Mean value
InitValueRandStCont1 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_student") Stddev value
InitValueRandStCont2 = ""
# {condition}(OperationType =="fill" and Fill_Operation=="rand_student") Df value
InitValueRandStCont3 = ""


#  {list_combo} (overwrite, append) Select writing mode
""" 
Different types of metadata writing mode
"""
ModeType = "overwrite"

# {hidden} Program name
"""This is the name of the program to be executed, dont change this!!!"""
ProgramName = "xmipp_metadata_utilities"

#{hidden} Usage of the program
Usage = """Perform several operations on metadata files."""

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------

from protocol_metadata_utilities import *
from protlib_xmipp import *

#        
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtMetadataUtils)
