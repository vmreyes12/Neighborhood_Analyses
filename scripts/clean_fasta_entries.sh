#!/bin/bash

# Remove empty entries | remove duplicates
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' $1 | awk '/^>/{f=!d[$1];d[$1]=1}f' > $2
