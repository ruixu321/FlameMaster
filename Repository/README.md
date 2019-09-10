## 1. What is FlameMaster?

FlameMaster is a C++ Program Package for 0D Combustion and 1D Laminar Flame 
Calculations.


## 2. Prerequisites

Operating system (OS) specific instructions for the prerequisites are in 
the directory `Repository/doc/`. 

For unix-like systems __`check your package manager`__ instead of installing
dependencies manually.  On Windows systems (except for
[cygwin](https://www.cygwin.com/)) flex and bison will be installed
automatically. Sundials can also be installed automatically (using the
`-DINSTALL_SUNDIALS=ON` option when you run cmake).

The following programs are required to install FlameMaster:

1. c and c++ compiler (C11/C++11 compliant) 
2. [cmake](http://www.cmake.org/download/)
3. [flex](http://sourceforge.net/projects/flex/files/)
4. [bison](https://www.gnu.org/software/bison/)
5. [sundials](https://computation.llnl.gov/casc/sundials/main.html)

## 3. Installation

### 3.1 Overview

The following is a brief introduction to CMake and an explanation of the steps 
required to install FlameMaster. You can skip this section if you already know 
CMake and build automation tools. Specific instructions for the installation 
based on the default configuration can be found in section 3.2.

FlameMaster is installed in three steps:

1. Configuration of the FlameMaster settings for your system using CMake
   (example command: `cmake ../Repository -DCMAKE_BUILD_TYPE=Release`).
2. Compilation using a build automation tool (example command: `make`).
3. Installation using a build automation tool (example command: `make install`). 
   Note that with the default configuration this step creates 
   sibling directories of `Repository/`.

- All prerequisites (for OS specific hints see `Repository/doc/`) need to be 
  found by CMake before you can build the FlameMaster project. If you installed
  the prerequisites into custom locations, CMake cannot know where to find the
  prerequisites and will complain. In that case, you must configure
  FlameMaster manually in the ccmake graphic user interface (from your build 
  directory run `ccmake ../Repository`). This interactive interface can also be 
  used to activate optional parts of the FlameMaster project or to change the 
  configuration (you can e.g. change the installation directory or include 
  debugging flags)

### 3.2 Installing FlameMaster with the default configuration

If not yet done, download and extract the FlameMaster package
or use `git clone` to get the repository. To use `git clone` either setup a 
[ssh-key](https://docs.gitlab.com/ee/gitlab-basics/create-your-ssh-keys.html) 
first and run

~~~~~~~~~~~{sh}
cd Your/Choice/FlameMaster
git clone git@git.rwth-aachen.de:ITV/FlameMaster.git Repository
~~~~~~~~~~~

or access the repository via `https` after you set your password

~~~~~~~~~~~{sh}
cd Your/Choice/FlameMaster
git clone https://git.rwth-aachen.de/ITV/FlameMaster.git Repository
~~~~~~~~~~~

Your username is either your github username or the first part of your
RWTH Aachen University email address (e.g. `max.mustermann`@rwth-aachen.de or 
`m.mustermann`@itv.rwth-aachen.de) depending on how you signed up.

From the command line execute the following (if 
[make](https://en.wikipedia.org/wiki/Make_(software)) is installed):

~~~~~~~~~~~{sh}
mkdir Build && cd Build
cmake ../Repository -DCMAKE_BUILD_TYPE=Release
make
make install
~~~~~~~~~~~

The command `cmake ../Repository -DCMAKE_BUILD_TYPE=Release` points CMake to 
the top level source code directory (`Repository/`) containing a CMakeLists.txt 
file and configures the FlameMaster project in the current directory (`Build/`). 
If the configuration finishes without errors, you can compile the FlameMaster 
package using `make`. This 'out-of-source build' leaves the `Repository/`
directory unchanged and allows you to set up multiple configurations. The 
command `make install` creates the installation directories of FlameMaster.

Note that CMake supports further 
[generators](https://cmake.org/cmake/help/v3.0/manual/cmake-generators.7.html), 
as for example [Microsoft Visual Studio](https://www.visualstudio.com/), if 
make is not available on your sytem.

### 3.3 Directory layout

Along with the introduction of CMake we decided to revise the directory layout 
of the project. This overview is supposed to help you to get started. After 
running the installation step there is the following default directory layout:

`FlameMaster/`
- `Repository/` permanent directory for the source code. __All changes__ of 
		the project files have to be done in this directory.
- `Build/` intermediate directory containing primarily __configuration__ 
		information ("CMake cache") and object files (*.o)
- `Run/` directory for running FlameMaster and its __examples__; Note that 
		files created during the installation will be overwritten if you 
		execute the installation again (`make install`). For permanent changes 
		modify instead the respective files in `Repository/examples/`
- `Bin/` directory containing __binaries__ and scripts that can be used to 
		set up your __environment__
- `Data/` directory that can be used as a default location to store binary 
		__input files__ like preprocessed mechanisms files (*.pre). This 
		directory is checked automatically by ScanMan and FlameMan via your 
		shell variable FM_DATA to find binary input files


## 4. Run FlameMaster

### 4.1. Environment

While not required to run FlameMaster it is more convenient to set up your 
environment first. Scripts to do this automatically are provided for bash, csh, 
and zsh:

~~~~~~~~~~~{sh}
source Your/Choice/FlameMaster/Bin/bin/Source.{yourshell}
~~~~~~~~~~~

You can make the setup permanent by adding the above line to your 
`~/.{yourshell}rc` file. The current setup can be verified and updated with the 
alias command `fmagain` which prints a list of all FlameMaster 
environment variables and alias commands. Note that manual changes in 
Source.{yourshell} will be overwritten if you rerun the installation. Instead, 
you should modify the file `Repository/config/SourceTemplate.{yourshell}.in` and 
run the installation command again (`make install`).

If you use none of the above shells it is also possible to run FlameMaster 
without setting up your environment or to configure your environment manually. 
A list of commands and variables can be found e.g. in `Bin/bin/Source.bash`.

### 4.2. Examples

Check out some examples in `FlameMaster/Run`. In each directory under Run, look 
at the README.md files, if present. These have descriptions of how to run the 
code and pre- and post-process the data. Also some additional helpful features 
are documented in the README.md files. 


## 5. Help

If you come across any problems during the installation, please send us the
complete output of the configuration step (`cmake ../Repository -DCMAKE_BUILD_TYPE=Release`) 
or the compilation step (`make`). 

Everyone is welcome to contribute to FlameMaster. You can for example report 
bugs by filing [issues](https://git.rwth-aachen.de/ITV/FlameMaster/issues). 
More comprehensive contributions require 
[elevated access](https://docs.gitlab.com/ee/user/permissions.html). You are 
welcome to contact us. Further information can be found in the
[contribution guide](./CONTRIBUTING.md).

