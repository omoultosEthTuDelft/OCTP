#!/bin/bash

outname=frdf_omp
# check if user specified $1 and it is a directory
if [ -z "$1" ]; then
    echo "Need at least one argument"
    echo "Usage: ./frdf_post.sh <input directory> <output file=frdf.dat>"
    exit 1
fi

if [ ! -d "$1" ]; then
    echo "Input directory $1 does not exist."
    exit 1
fi

# if $2 exists, make it the output name
if [ ! -z "$2" ]; then
    outname = $2
fi

# check if $2 exists in $1, if it does, abort
if [ -f "$1/$outname" ]; then
    echo "Output file $1$outname already exists."
    exit 1
fi


cd $1
../../rdf_post/frdf_omp dumpfile temperature.dat $outname 1000 1000
