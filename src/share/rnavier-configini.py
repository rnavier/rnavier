#!/usr/bin/env python

import ConfigParser
import shutil
import re
import sys
import os
import argparse

# Opens a configuration file for reading
def find_in_path(name, path):
    """Search PATH for a binary.

    Args:
      name: the filename to search for
      path: the path  ['./', './path/to/stuff']

    Returns:
      The abspath to the fie or None if not found.
    """
    for dir in path:
        binpath = os.path.join(dir, name) 
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None

# Recursively searches for include filenames, keeping track of order
def find_include_filenames(filenames, files, path) :
    try:
        filein = open(filenames[-1], "r") 
    except IOError:
        print "Error: Could not open ", filenames[-1], " requested in INCLUDE directive"
        sys.exit()
    files.append(filein)
    pattern = re.compile('##INCLUDE\s+([\w.\/]+)') 

    lines = files[-1].readlines()
    lines.reverse()
    for line in lines:
        m = pattern.match(line) 
        if m :
            name = find_in_path(m.group(1), path)
            if not name:
                print "Requested filename", m.group(1), "is not in path."
                sys.exit()

            if  not name in filenames:
                filenames.append(name)
                find_include_filenames(filenames, files, path)
                continue
            else:
                print "Warning! Multiple INCLUDE directives are adressing the same file", m.group(1)
    files.pop().close()

#
# Parse arguments
#
parser = argparse.ArgumentParser(description = "configureini.py reads filename.cfg and produces filename.ini, which can be used as input for the hydrodynamics code")

parser.add_argument('-I', '--include', metavar='path/to/include', action='append', help='include path for configuration files')
parser.add_argument('filename', metavar='filename.cfg', help='the input configuration file')
args = parser.parse_args() 


# Construct the search path [./, included paths, 'RNAVIERDATA/defaults']
path = []
path.append(os.getcwd())
if args.include:
    for name in args.include:
        if os.path.exists(os.path.abspath(name)) :
            path.append(os.path.abspath(name))

if 'RNAVIERDATA' in os.environ:
    path.append(os.path.join(os.environ['RNAVIERDATA'], 'defaults'))
else:
    print "Warning: environment variable RNAVIERDATA is not defined."


## Checkthat filename has a cfg extension and is in path
name = find_in_path(args.filename, path)
if not name:
    print "Could not find configuration file", args.filename,"in search path!" 
    sys.exit()

(root, ext) = os.path.splitext(name)
if (ext != '.cfg') :
    print "Error: input filename must have a .cfg extension, i.e. filename.cfg"
    sys.exit()

## Use find_include_filenames to determine which files are actually included
include_filenames = [name]
find_include_filenames(include_filenames, [], path)
include_filenames.reverse()
#print include_filenames


# Construct a configuration parser based on these files
config = ConfigParser.ConfigParser()
config.optionxform = str
config.read(include_filenames) 

## Write filename.ini
name = find_in_path(args.filename,path)
(root, ext) = os.path.splitext(name)
filename_ini = root + '.ini'

fileout = open(filename_ini, "w")
print  >> fileout, "# Automatically generated from:"
for name in include_filenames:
    print >>fileout, "#", name

# Write the ini file
config.write(fileout)

print "Successfully created", filename_ini
