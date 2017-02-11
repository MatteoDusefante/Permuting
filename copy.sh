#!/bin/bash


if [ $1 = "-fc" ]; then
    scp -P 41122 madu@sssgate.itu.dk:~/output.txt  $2
    exit 0
fi

if [ $1 = "-tc" ]; then
    for i in $( ls ); do
        scp -P 41122 $i madu@sssgate.itu.dk:~/
    done
    exit 0
fi

if [ $1 = "test" ]; then
     scp -P 41122 papi_test.cpp papils.h madu@sssgate.itu.dk:~/
     exit 0
fi


scp -P 41122 papils.h Makefile papi_omp.cpp madu@sssgate.itu.dk:~/

if [ $# -eq 0 ]; then
    exit 0
fi

if [ $1 = "spmxv" ]; then
     scp -P 41122 spmxv.cpp utl.h madu@sssgate.itu.dk:~/
     exit 0
fi

if [ $1 = "bench" ]; then
    scp -P 41122 benchmark.cpp benchmark.h best_bucket.cpp utils.h algorithms.h madu@sssgate.itu.dk:~/
    exit 0
fi

if [ $1 = "cuda" ]; then
    scp -P 41122 cuda.cu cutils.cuh utils.h madu@sssgate.itu.dk:~/
    exit 0
fi

if [ $1 = "perm" ]; then
    scp -P 41122 permutation.cpp permutation.h utils.h madu@sssgate.itu.dk:~/
    exit 0
fi
