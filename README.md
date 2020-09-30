# rlibm-generator
This tool is used to generate the [rlibm](https://github.com/rutgers-apl/rlibm) math library. For more information behind the techniques in **rlibm-generator**, please refer to this [paper](https://arxiv.org/pdf/2007.05344.pdf).

## Installation

### Prerequisite:
* GMP : Needed by Soplex 4.0.1
* [Soplex 4.0.1](https://soplex.zib.de/index.php#download) : Please read below on how to install Soplex 4.0.1 to use with rlibm-generator
* MPFR
* [SoftPosit](https://gitlab.com/cerlane/SoftPosit)

### Installation instruction:
1. Make sure you have GMP and MPFR installed.

2. Install SoftPosit according to the instructions from the [SoftPosit website](https://gitlab.com/cerlane/SoftPosit)

3. Install [Soplex 4.0.1](https://soplex.zib.de/index.php#download)

  1. Download Soplex version 4.0.1 source code.
  
  2. Untar Soplex
  ```
  tar -xvf soplex-4.0.1.tar
  ```
  
  3. cd into soplex directory and create a build directory inside
  ```
  cd soplex-4.0.1
  mkdir build
  ```
  
  4. cd into build directory and build soplex. Make sure you have GMP installed prior to this step.
  ```
  cd build
  cmake ..
  make
  ```
  
4. Create 2 environment variables to specify where softposit and soplex is installed:
```
export SOFTPOSITPATH=<path_to_softposit_directory>
export SOPLEXPATH=<path_to_soplex-4.0.1_directory>
```

5. To run the math library function generators, use the script
```
cd rlibm-generator
./runFunctioonGen.sh
```

## More Detail

If you would like to generate the functions individually, the function generation scripts are in functiongen directory. Functiongen directory is organized into each representation (bfloat16, float, posit16, etc) and each directory is further organized into the method of how we generate the function. There are two ways we generate the functions:
1. Full: We use ALL intervals when generating the LP query and feed it to the LP solver to generate the polynomial
2. Sampling: We sample some intervals, create LP query based on the sample.

Inside of full or sampling directory, there are cpp files that generates the polynomial for each elementary function and a makefile. You can use the makefile to compile the program and run each program that generate the polynomial for a single elementary function.
