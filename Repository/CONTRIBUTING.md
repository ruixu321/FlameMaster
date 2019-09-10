# How to contribute

Everyone is welcome to contribute to `FlameMaster`. 
Contributing doesn’t just mean submitting merge requests — 
there are many different ways for you to get involved, 
including writing tutorials in the [Wiki](https://git.rwth-aachen.de/ITV/FlameMaster/wikis/home), 
contribute new or improve existing 
[tests cases](https://git.rwth-aachen.de/ITV/FlameMaster/pipelines), reporting 
bugs by creating [issues](https://git.rwth-aachen.de/ITV/FlameMaster/issues), 
and participating in the `FlameMaster` evolution process.

There are a few guidelines that we need contributors to follow so that we can 
have a chance of keeping on top of things.

# Reports bugs

Reporting bugs is a great way for anyone to help improve `FlameMaster`.
Submit a new [issue](https://git.rwth-aachen.de/ITV/FlameMaster/issues)
for the bug, assuming one does not already exist.
When filing a new `FlameMaster` bug, please include the following:

- A concise description of the problem.

- A reproducible test case. If the issue occurs when you run the 
`FlameMaster` executables, paste the call into the issue’s description field. 
Also provide alias commands and environment varialbes that you used (e.g. the 
output of the `fmagain` command). 

- Attach all required files to run your test case, if they are not 
already a part of the repository. These can be:
    - Input file (`FlameMaster.input`)
    - Start profile
    - Mechanism file
    - Theromdynamic data file
    - Transport data file

- If necessary, attach the redirected screen output (as `screen.log`) and 
generated output files.

- A description of the environment that reproduces the problem. If possible, 
include the git hash (`git rev-parse --verify HEAD`), as well as information
about your OS and platform.

- Details about the installation. You can for example provide the commands
run to install `FlameMaster`.

# Comprehensive Contributions

Comprehensive contributions such as merge requests or 
[wiki articles](https://git.rwth-aachen.de/ITV/FlameMaster/wikis/home) require 
[elevated access](https://docs.gitlab.com/ee/user/permissions.html). Please 
contact us for further information.