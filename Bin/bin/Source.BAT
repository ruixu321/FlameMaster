:: List of alias commands (some of the command are redundant):
::	- FlameMan 		(to execute FlameMan)
::	- FlameMaster 	(to execute FlameMan)
::	- ScanMan 		(to execute ScanMan)
::	- CreateBinFile (to execute CreateBinFile)
::	- ListTool 		(to execute ListTool including the -M option)
::	- LT 			(to execute ListTool including the -M option)

@ECHO off
IF "%1"=="" (
@ECHO. 
@ECHO    Setting up doskeys to run FlameMaster
@ECHO.
)
@SET DA_ROOT="%~dp0"
:: @ECHO %DA_ROOT%
:: remove double quotes
@SET DA_ROOT=%DA_ROOT:"=%
IF "%1"=="" (@ECHO    Script location is "%DA_ROOT%")
::@DOSKEY DR_ROOT_BIN=Bin/bin
::@DOSKEY DR_ROOT_LIB=Bin/lib
::@DOSKEY DA_ROOT_DATA=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/../Data

@SET DR_ROOT_BIN=Bin/bin
@SET DR_ROOT_LIB=Bin/lib
@SET DR_SRC_SCRIPT_LOCATION=../..\
@SET DA_ROOT_DATA=/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/../Data
IF EXIST "%DA_ROOT_DATA%\." (
	IF "%1"=="" (@ECHO        Data folder "%DA_ROOT_DATA%" exists)
) ELSE (
	@SET DA_ROOT_DATA="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%Data"
)

@SET DA_ROOT_DATA=%DA_ROOT_DATA:"=%

IF "%1"=="" (@ECHO        Checking possible location if you installed the binaries)
IF EXIST "%DA_ROOT_DATA%\." (
	IF "%1"=="" (@ECHO        Data folder "%DA_ROOT_DATA%" exists)
) ELSE (
	IF "%1"=="" (@ECHO        Creating FlameMaster data directory "%DA_ROOT_DATA%")
	@MKDIR "%DA_ROOT_DATA%"
)

@SET DR_ROOT_RUN=Run

@SET FM_BIN=%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%

IF "%1"=="" (
	@ECHO    Creating doskeys...
	@ECHO        FlameMan.exe -^> FlameMan
)
@DOSKEY FlameMan="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\FlameMan.exe" $*
@DOSKEY FlameMaster="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\FlameMan.exe" $*
IF "%1"=="" (@ECHO        ScanMan.exe -^> ScanMan)
@DOSKEY ScanMan="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\ScanMan.exe" $*
IF "%1"=="" (@ECHO        CreateBinFile.exe -^> CreateBinFile)
@DOSKEY CreateBinFile="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\CreateBinFile.exe" $*
IF "%1"=="" (@ECHO        ListTool.exe -^> ListTool)
@DOSKEY ListTool="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\ListTool.exe" -M $* 
@DOSKEY LT="%DA_ROOT%%DR_SRC_SCRIPT_LOCATION%%DR_ROOT_BIN%\ListTool.exe" -M $*
IF "%1"=="" (
	@ECHO.
	@ECHO    Setting FlameMaster data path...
)

@SET FM_DATA="%DA_ROOT_DATA%"
@SET FM_DATA=%FM_DATA:"=%
IF "%1"=="" (@ECHO    FM_DATA="%DA_ROOT_DATA%")


:: Set up some more linux-like behavior for windows!
::@DOSKEY ls=dir $*
::@DOSKEY cd=cd $1$Tdir
::@DOSKEY clear=cls
::@DOSKEY rm=del $*
::@DOSKEY cp=copy $*

IF "%1"=="" (
	@CD "%DR_SRC_SCRIPT_LOCATION%"
	@ECHO.
	@ECHO    done!
	@ECHO.
	@ECHO    Enjoy using FlameMaster 4.0.0
	@ECHO. 

	@cmd /K
)
