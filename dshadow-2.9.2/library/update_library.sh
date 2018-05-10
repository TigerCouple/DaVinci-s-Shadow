#!/bin/bash

#get the most recent library copy
svn update

#get the log for the library
svn log -v | sed -e 's/^[ \t]*//' | sed -e 's/^[+]/*/' | head -n -1> change_log_dvrc.txt

#get the log for the core
svn log -v http://oss.mars.asu.edu/svn/davinci/davinci/trunk | sed -e 's/^[ \t]*//' | sed -e 's/^[+]/*/' | head -n -1> change_log_core.txt
svn log -v http://oss.mars.asu.edu/svn/davinci/davinci/trunk/version.h | sed -e 's/^[ \t]*//' | sed -e 's/^[+]/*/' | head -n -1> change_log_version.txt

#make the update_library.txt download file
find -L . -name "*.*" | grep -v "~" | grep -v "#" | grep -v "svn" | cut -c 3- | tail -n +2 | grep -v "update_library" > update_library.txt
