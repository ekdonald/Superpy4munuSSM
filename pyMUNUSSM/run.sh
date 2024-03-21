#!/bin/bash

LD_LIBRARY_PATH=../MultiNest-master/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
#python3 SuperPy.py

python3 packageinfo.py
mpiexec -np 1 python3 SuperPy.py
