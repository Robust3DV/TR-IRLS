# TR-IRLS: Robust Fitting with Truncated Least Squares
Source code of the [3DV'21](https://3dv2021.surrey.ac.uk/) paper:

[Robust Fitting with Truncated Least Squares: A bi-level Optimization Approach.](https://www.researchgate.net/publication/356834866_Robust_Fitting_with_Truncated_Least_Squares_A_Bilevel_Optimization_Approach)


## Obtaining the source code:

```
git clone https://github.com/Robust3DV/TR-IRLS.git 
```


## Compilation:
The source code is written in C++ based on the [SSBA](https://github.com/chzach/SSBA/tree/master/SSBA-4.0) library. It has been tested on an Ubuntu machine. 

In order to compile the source code, from the repository's main folder, execute the following commands:


``` 
mkdir build
cd build
cmake ..
make -j8 

```

## Usage:
If the compilation is successful, an executable file (***bundle_trunc***) will be created. The syntax is

```
./bundle_trunc <path_to_bal>
```
where **path_to_bal** is the path to the an instance of bundle adjustment (in [BAL's format](https://grail.cs.washington.edu/projects/bal/)).

An example input file is provided in the Datasets folder. To run the example, execute the following command from the ***build*** folder:
```
./bundle_trunc ../Datasets/problem-16-22106-pre.txt
```


