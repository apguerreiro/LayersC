LayersC
=====

This software implements the LayersC algorithm for the Hypervolume Subset Selection Problem (HSSP) in three dimensions.

**NOTE**
This is a beta version of the implementation of the algorithm proposed in [1]. The code is being prepared for release. The code depends on a modified versions of HVC and EAF3D algorithms which will be added/updated soon.

License
--------


Except where indicated otherwise in individual source files, this software is Copyright © 2020-2022 Andreia P. Guerreiro.

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Appropriate reference to this software should be made when describing research in which it played a substantive role, so that it may be replicated and verified by others. The algorithm which this software implements is described in detail in [1]. 


Dependencies
--------
- [HVC](https://github.com/apguerreiro/HVC)
- [EAF](https://eden.dei.uc.pt/~cmfonsec/aft-0.95.tar.gz)
- [gHSS](https://github.com/apguerreiro/gHSS)
- CPLEX

Building
--------



In GNU/Linux, the program can be compiled from source by invoking:

    make


General Usage
-------------


**SYNOPSIS** 

    layers [OPTIONS] [FILE...]
    
        
**DESCRIPTION**

Compute a solution to the Hypervolume Subset Selection problem for the data set(s) in FILE(s).

With no FILE, read from the standard input.

**COMMAND LINE OPTIONS**

			Options:
			 -h, --help          print this summary and exit.                          
			     --version       print version number and exit.                        
			 -v, --verbose       print some information (time, coordinate-wise maximum 
					     and minimum, etc)                                     
			 -q, --quiet         print just the hypervolume (as opposed to --verbose). 
			 -u, --union         treat all input sets within a FILE as a single set.   
			 -r, --reference=POINT use POINT as reference point. POINT must be within  
					     quotes, e.g., "10 10 10". If no reference point is  
					     given, it is taken as the coordinate-wise maximum of  
					     all input points.                                     
			 -k, --subsetsize=k  select k points (a value between 1 and n, where n is  
					     the size of the input data set. The default is n/2)   
			 -f, --format=(0|1) output format                                          
					     (0: print selected points and the hypervolume indicator)
					     (1: print selected points)                            
			 -g, --greedystart   start with a greedy solution         
                               	


Usage
-----

**Run**

The program reads sets of points in the file(s) specified in the command line:

    ./layers data


In input files, each point is given in a separate line, and each coordinate within a line is separated by whitespace. An empty line, or a line beginning with a  hash sign (#), denotes a separate set.


Sets in an input file may be treated as a single set by using option `-u`: 

    ./layers data -u


A reference point can be set by giving option `-r`.

    ./layers -r "10 10 10 10" data

 If no reference point is given, the default is the coordinate-wise maximum of all input points in all files.

For the other options available, check the output of `./layers --help`.
 
By default, the default subset size *k* is set to *n/2* where *n* is the data set size. The subset size can be explicitly specified using option `-k`:
 
    ./hvc -r "10 10 10 10" data -k 10
    
The use of the option `-g` is recommended, for better performance as the algorithm starts with an approximated solution.

For the remainder options available, check the output of `./hvc --help`.



References
----------



[1] Guerreiro, Andreia P.; Manquinho, Vasco; Figueira, José Rui. "Exact hypervolume subset selection through incremental computations". Computers & Operations Research 136 (2021): 105471. http://dx.doi.org/10.1016/j.cor.2021.105471.




