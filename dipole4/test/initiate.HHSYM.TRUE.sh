#!/bin/bash

rm -f *.dat

echo "compute energies (HHSYM=TRUE)"
../main.exe.HHSYM.TRUE > output.txt

rm -f h6a.dat h6s.dat h6.dat

echo "transform to primite basis"
../transform2PB/block/main.exe > out.transf.txt

echo "done"
