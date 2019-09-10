# Prerequisites: MacOS X

You should be able to install all prerequisites using 
[homebrew] (http://brew.sh) or [macports] (https://www.macports.org). The 
package `gcc` (containing `gfortran`) is optional. Other prerequisites are 
included in the Mac OS X command line tools that are required to install 
`homebrew`.

## Homebrew (03/27/2017)

With homebrew the following command installs all prerequisites:
~~~~~~~~~~~{sh}
brew install git cmake homebrew/science/sundials gcc

or 

brew install git cmake homebrew/science/sundials gcc // Added by Rui 10/26/2018
~~~~~~~~~~~

## Macports (03/27/2017)

~~~~~~~~~~~{sh}
sudo port install cmake git sundials gcc49
~~~~~~~~~~~
