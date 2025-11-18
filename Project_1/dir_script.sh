#!/bin/bash

#take path as command line input
command_line="$1"
#make the specified directory the working
mkdir -p "$command_line"
#enter the directory
cd "$command_line"
#print the abs path
pwd


#./dir_script.sh "/Users/ark/Documents/GT/BINF/Fall/Program Bioinformatics"