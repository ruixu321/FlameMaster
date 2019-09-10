# Prerequisites: Windows

On Windows systems there are generally two options to install FlameMaster. 
The first option is to install FlameMaster natively using for e.g. 
[Microsoft Visual Studio](https://www.visualstudio.com/)
or [cygwin](https://www.cygwin.com/). The second option is to install FlameMaster 
in a [containerized (docker)](https://en.wikipedia.org/wiki/Docker_(software)) or 
[virtualized (e.g. Virtual Box)](https://en.wikipedia.org/wiki/Virtualization) system. 
The advantage of the first option is reduced overhead, while the advantage of the 
second option is, that FlameMaster is well validated for unix-like environments.
Using the first option the installation procedure can be slightly different from 
the default installation as `make` is not the most common build tool on Windows. 
If you choose the second option, the procedure to build and install FlameMaster 
is the same as on unix-like systems. 

## Native builds on Windows

If you do not use [cygwin](https://www.cygwin.com/), flex and bison 
will be installed automatically by FlameMaster. Otherwise, make sure that
the [cygwin](https://www.cygwin.com/) installation includes these programs.
In both cases, it is recommended to let FlameMaster 
build sundials from source. The sundials installation requires
[python](https://www.python.org/downloads/) to be installed. Both commands 
`python` and `cmake` must be available from the command line which means the
locations of the executables must be added to your `PATH` environment variable. 
Note that this might not be the default option and must be selected when you 
install the programs. You have to reboot to update the `PATH` variable. 

Run the following commands to install the FlameMaster package:
~~~~~~~~~~~{sh}
mkdir Build && cd Build
cmake ../Repository -DCMAKE_BUILD_TYPE=Release -DINSTALL_SUNDIALS=ON
# only with MinGW...
mingw32-make
mingw32-make install
~~~~~~~~~~~

With [cygwin](https://www.cygwin.com/) or [MinGW](http://www.mingw.org/)
you will probably use `make` or `mingw32-make` to build and install FlameMaster. 
It is also possible to build FlameMaster in an 
[IDE](https://en.wikipedia.org/wiki/Integrated_development_environment) 
(e.g. [Eclipse](https://eclipse.org/ide/) 
or [Microsoft Visual Studio](https://www.visualstudio.com/)). The specific 
steps are generally similiar to the ones decribed here for
[Microsoft Visual Studio](https://www.visualstudio.com/).
With [Microsoft Visual Studio](https://www.visualstudio.com/) you should 
set the configuration to `Release` and build the `INSTALL` solution. 
After the installation in [Microsoft Visual Studio](https://www.visualstudio.com/), 
a simple batch file is available in the directory `Bin/bin` that allows you to
run the FlameMaster executables via `doskeys`.

## None native installation (with Docker)

This is a brief overview of what you need to know to run FlameMaster on a unix-like system as 
a Windows user. [Virtualizations (e.g. Virtual Box)](https://www.virtualbox.org/) will not be 
discussed.

Please look up any Linux or Unix tutorial to familiarize yourself with the linux
command line syntax. Many commands are similar to Windows, but some are 
different (`dir` = `ls`, `copy` = `cp`, `ren` = `mv`, etc.). The shorthand for your 
home directory is `~`. Directory names are separated using `/` as opposed to `\`.

To install FlameMaster on Windows computers, you can use [Docker](https://www.Docker.com/).
With Docker you can create lightweight images and run them in containers, similar to a 
[virtual machine](https://en.wikipedia.org/wiki/Virtual_machine). We use such an 
image (Linux, Ubuntu based) to provide all prerequisites and to run FlameMaster. 
You require at least 500 MB of space on your hard drive for the FlameMaster
image. See documentation of Docker for any Docker installation and usage issues.

To set up the FlameMaster-docker-container execute following steps:
1. [Download and install Docker](https://www.Docker.com/products/Docker#/windows)
2. If not yet done, download and extract the FlameMaster package (to `Your/Choice/FlameMaster`). 
Run your favorite command line and navigate to `Your/Choice/FlameMaster/Repository`
3. Make sure docker is running and execute:

~~~~~~~~~~~{sh}
docker build --pull -f docker/ubuntu.df -t ubuntufm .
docker run --name fmrun -it ubuntufm bash
~~~~~~~~~~~
	
The command `docker build` creates an image from the Dockerfile `docker/ubuntu.df`. 
This image includes the FlameMaster package preinstalled in the directory `/FlameMaster`. 
The option `-t` is used to tag the image as `ubuntufm`.
The `docker run` command creates a new container (or "instance") of the flamemaster image 
called `fmrun`. The `-it` flags start an interactive bash shell in the container. 

You can now use the FlameMaster package in the same way as on Ubuntu. The container stops 
running if you exit the bash shell. You should use the `docker run` command only once to 
create the container. From the second time on, restart your container by typing the following
in the command line

~~~~~~~~~~~{sh}
docker start -i fmrun
~~~~~~~~~~~

For additional information visit the 
[Docker documentation](https://docs.Docker.com/engine/reference/commandline/cli/) 
or type `docker --help`. For further information about the installtion you might 
want to review the section "Installtion" in the top level README.md file.
