#!/bin/bash

rm -f *.dat

echo "compute energies (HHSYM=FALSE)"
../main.exe.HHSYM.FALSE > output.txt

rm -f h6a.dat h6s.dat h6.dat

echo "transform to primite basis"
../transform2PB/block/main.exe > out.transf.txt

echo "done"
