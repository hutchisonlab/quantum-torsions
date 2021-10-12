#!/bin/bash

bzip2 -dk *.bz2

for i in *.out
do
    mv -v "${i}" "${i%.*}_dft.${i##*.}"
done

obabel *_dft.out -osdf -m

