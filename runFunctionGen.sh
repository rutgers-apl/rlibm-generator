#!/bin/bash

echo -e "\033[1mCompiling bfloat16 function generators...\033[0m"
cd functiongen/bfloat16/full
make -s clean
make -s
echo -e "\033[1mGenerating bfloat16 functions... \033[0m"
./runAll.sh
make -s clean
cd ../../..

echo -e "\033[1mCompiling posit16 function generators... \033[0m"
cd functiongen/posit16/full
make -s clean
make -s
echo -e "\033[1mSynthesizing posit16 functions \033[0m"
./runAll.sh
make -s clean
cd ../../..

echo -e "\033[1mCompiling float function generators... \033[0m"
cd functiongen/float/sampling
make -s clean
make -s
echo -e "\033[1mSynthesizing float functions \033[0m"
./runAll.sh
make -s clean
cd ../../..
