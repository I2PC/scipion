#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
""" 
Describe your run here...
"""
#-----------------------------------------------------------------------------
# {section} Run 
#-----------------------------------------------------------------------------
# Run name:
""" 
This will identify your protocol run. It need to be unique for each protocol. 
You could have <run1>, <run2> for protocol X, but not two run1 for same protocol. 
This name together with the protocol output folder will determine the working
directory for this run.
"""
RunName = "run_001"

# {list}(Resume, Restart) Run behavior
""" 
Resume from the last step, restart the whole process 
or continue at a given step or iteration.
"""
Behavior = "Resume"