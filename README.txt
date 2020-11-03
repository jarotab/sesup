=================================================
SESUP
Testing tool for state estimation algorithms for systems with uncertain parameters.
Author: Jaroslav Tabacek (tabacjar@fel.cvut.cz)

=================================================
ABOUT
-------------------------------------------------
The SESUP application is created for testing and comparing state estimation alogrithms for systems with uncertain parameters.
This version supports the results of paper 
J. Tabacek and V. Havlena, "Reduction of stochastic prediction error sensitivity to parameters in Kalman filter," Journal of Franklin Institute, submitted, Nov. 2020

=================================================
FILES
-------------------------------------------------
1. SesupInstallerWeb.exe 	 - standalone app installer for Windows
1. SesupInstallerWeb.install - standalone app installer for Linux
2. SESUP.mlappinstall    	 - MATLAB app installer for MATLAB (R2019b and later)
3. data.zip				 	 - saved simulation data which are presented in the paper (can be imported into the app)
4. readme.txt            

=================================================
INSTALLATION
-------------------------------------------------
There are two ways how to use the app:

1. Standalone application for Windows/Linux
	- MATLAB license and installation are NOT required.
	- WINDOWS (Tested in Windows 10 Enterprise)
		- install app by running 'SesupInstallerWeb.exe' and follow the instructions
		- SESUP application will be installed together with MATLAB Runtime	
		- start the app by locating and running 'SESUP.exe'		
	- LINUX (Tested in Ubuntu 20.04) 
		- install app by running "SesupInstallerWeb.install" and folllow the instructions (don't forget to configure environment variables)
		- SESUP application will be installed together with MATLAB Runtime		 
		- start the app by locating and running "run_SESUP.sh"		

2. MATLAB App
	- install the app directly to your MATLAB (R2019b or later) or to MATLAB Online
	- Windows and Linux computers are supported
	- run 'SESUP.mlappinstall' to install the app to the MATLAB
	- the app will be added to MATLAB App ribbon after installation
	- start the app by clicking to the application icon at MATLAB app ribbon  

=================================================
LIMITS
-------------------------------------------------
1. Standalone application for Windows/Linux (limited functionality)
	A. Tested in Windows 10 Enterprise / Ubuntu 20.04
	B. Adding your custom model for testing is not available
	C. Editing and simulating will show the warning that code generation failed. This message can be ignored, since the pregenerated code for the default system described in the paper is included.

2. MATLAB App 
	A. Tested in MATLAB R2019b, R2020a and R2020b (not recommended). There are regular behavior changes in released MATLAB App Designer, therefore, the app might not work in newer versions.
	B. Importing and browsing data does not require any additional toolbox mentioned below.	
	C. Requirements for running simulations:
		- Casadi toolbox in MATLAB path. Tested with version 3.0 and 3.5.5.
		- MATLAB supported C compiler. 
			Windows: Install via support package "MATLAB Support for MinGW-w64 C/C++ Compiler"
			Linux: GCC C/C++ usually comes with distribution. If not, it can be installed by command "sudo apt-get install gcc"
		- (Optional) "Parallel Computing Toolbox" can significantly speed up simulations. It is not required, but strongly recommended when running more Monte Carlo simulation.

3. MATLAB Online
	- The limitations of standalone app are combined with the limitations of MATLAB app
		- Unix-based machine - see 2C: linux version of Casadi toolbox must be in path
    	- Code compilation is not available - see 1B and 1C 