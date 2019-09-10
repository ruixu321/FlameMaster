# Prerequisites: Linux

## Ubuntu (03/27/2017)

You get everything you need with aptitude although the provided cmake of the 
`14.04 Ubuntu LTS version` is ancient. When installing sundials you must choose
 the sundials __dev__ packages. The regular sundials package doesn't create symbolic links 
in `/usr/lib` which means cmake cannot find the libraries. The command to 
install all prerequisites is:

~~~~~~~~~~~{sh}
sudo apt-get install git cmake flex bison libsundials-serial-dev g++ gfortran cmake-curses-gui
# cmake ...
~~~~~~~~~~~

`gfortran cmake-curses-gui` are optional. 

## Cent OS (03/27/2017)

~~~~~~~~~~~{sh}
yum install -y \
	bison \
	clang \
	clang++ \
	flex \
	gcc-c++ \
	gcc-gfortran \
	git \
	make \
	cmake \
	python

# cmake ...
~~~~~~~~~~~

## openSUSE (03/27/2017)

The recommended installation steps are:
~~~~~~~~~~~{sh}
zypper --no-gpg-checks refresh && zypper --non-interactive install \
	bison \
	cmake \
	flex \
	gcc \
	gcc-c++ \
	gcc-fortran \
	git 
cmake ../Repository -DCMAKE_BUILD_TYPE=Release -DINSTALL_SUNDIALS=ON 
make -j4 install
~~~~~~~~~~~ 

The above steps build sundials from source. If there are any 
issues with the sundials installation you can also try the 
sundials package from `zypper`.

The `openSUSE Leap` provided sundials installation does not inlcude `IDAS`. 
When installing sundials you must choose the sundials __dev__ packages. 

An installation that compiles exclusively the FlameMaster package from source looks like this:
~~~~~~~~~~~{sh}
zypper --non-interactive clean && zypper --no-gpg-checks --non-interactive ref && zypper --non-interactive up --no-confirm
# Choose the correct repository for your opensuse version
# https://software.opensuse.org/download.html?project=science&package=sundials
zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_Leap_42.2/science.repo
zypper --no-gpg-checks refresh && zypper --non-interactive install \
	bison \
	cmake \
	flex \
	gcc \
	gcc-c++ \
	gcc-fortran \
	git \
	sundials-devel 

cmake ../Repository -DCMAKE_BUILD_TYPE=Release -DSUNDIALS_NO_IDAS_SEARCH=ON 
make -j4 install
~~~~~~~~~~~

The option `SUNDIALS_NO_IDAS_SEARCH` was specifically implemented for this purpose and is 
generally not recommened.
