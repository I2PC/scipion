# -*- coding: utf-8 -*-

# the URLDownload-Builder can be download any data from an URL into a target file
# and can replace the target file name with the URL filename (the setting variable
# within the environment object is a boolean type with the name "URLDOWNLOAD_USEURLFILENAM", 
# default setting replaces the target name with the URL filename)

# JMRT 11-Nov-2013
# Added the Var URLDOWNLOAD_DIRECTORY to specify the output download dir

import os
import urllib2, urlparse
import SCons.Builder, SCons.Node, SCons.Errors


# define an own node, for checking the data behind the URL,
# we must download only than, if the data is changed, the
# node derivates from the Python.Value node
class URLNode(SCons.Node.Python.Value) :

    # overload the get_csig (copy the source from the
    # Python.Value node and append the data of the URL header
    def get_csig(self, calc=None): 
        try: 
            return self.ninfo.csig 
        except AttributeError: 
            pass 
        
        try :
            response = urllib2.urlopen( str(self.value) ).info()
        except Exception, e :
            raise SCons.Errors.StopError( e )
            
        contents = ""
        if "Last-Modified" in response :
            contents = contents + response["Last-Modified"]
        if "Content-Length" in response :
            contents = contents + response["Content-Length"]
        if not contents :
            contents = self.get_contents() 
        self.get_ninfo().csig = contents 
        return contents 



# creates the downloading output message
# @param s original message
# @param target target name
# @param source source name
# @param env environment object
def __message( s, target, source, env ) : 
    print "downloading [%s] to [%s] ..." % (source[0], target[0])


# the download function, which reads the data from the URL
# and writes it down to the file
# @param target target file on the local drive
# @param source URL for download
# @@param env environment object
def __action( target, source, env ) :
    try :
        output_file = str(target[0])
        stream = urllib2.urlopen(str(source[0]))
        file   = open( output_file, "wb" )
        file.write(stream.read())
        file.close()
        stream.close()
    except Exception, e :
        raise SCons.Errors.StopError( e )


# defines the emitter of the builder
# @param target target file on the local drive
# @param source URL for download
# @param env environment object
def __emitter( target, source, env ) :
    # we need a temporary file, because the dependency graph
    # of Scons need a physical existing file - so we prepare it
    target[0].prepare()
    
#    output_dir = env.get("URLDOWNLOAD_DIRECTORY", '.')

    if env.get("URLDOWNLOAD_USEURLFILENAME", False) :
        try :
            url = urlparse.urlparse( urllib2.urlopen( str(source[0]) ).geturl() )
        except Exception, e :
            raise SCons.Errors.StopError( e )
        target = [url.path.split("/")[-1]]
    
    output_dir = env.get("URLDOWNLOAD_DIRECTORY", '.')
    output_file = os.path.join(output_dir, str(target[0]))
    
    return output_file, source




# generate function, that adds the builder to the environment,
# the value "DOWNLOAD_USEFILENAME" replaces the target name with
# the filename of the URL
# @param env environment object
def generate( env ) :
    env["BUILDERS"]["URLDownload"] = SCons.Builder.Builder( action = __action,  emitter = __emitter,  target_factory = SCons.Node.FS.File,  source_factory = URLNode,  single_source = True,  PRINT_CMD_LINE_FUNC = __message )
    env.Replace(URLDOWNLOAD_USEURLFILENAME =  True )

# existing function of the builder
# @param env environment object
# @return true
def exists(env) :
    return 1
