### Created by wffx 

	Email : shuori2014@gmail.com

### base.hpp

This file is the falsework for Ring-LWE, including the class for polynomial calculation, Gaussian sampling and so on.
And it can be applied into other implemention based on Ring-LWE. 


### Dependence:
	
	1. NTL

	Follow the tour of the ntl website 
	Make true gmp and gf2x installed before install ntl.
	https://www.shoup.net/ntl/doc/tour-unix.html

    ./configure DEF_PREFIX=/usr/local PREFIX=$(DEF_PREFIX) SHARED=on NTL_GF2X_LIB=on 


### Complie command:
g++ -g -O2 -std=c++11 -pthread -march=native demo.cpp -o demo -lntl -lgmp -lm

or

make

### Note for my implementation:

1. There is still something wrong in the bound of parameters. to be continued ... 

