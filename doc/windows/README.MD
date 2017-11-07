# Windows
## Stuff To Install
Getting ZPIC up under windows is just a matter of a few easy steps. First, install the _Visual Studio_ - the standard Windows Programming environment from Microsoft - if you don't already have it. There is full featured free version at 
```
https://www.visualstudio.com/vs/community/
```
Don't worry about the installation options.. just choose the default ones and wait (it may take a bit of time) for the installation to finish. 

You need to install one other piece of software - _Cmake_. _CMake_ is a cross-platform building tool that makes building ZPIC easier on Windows. You can get it from:
```
https://cmake.org/files/v3.9/cmake-3.9.4-win64-x64.msi
```
Again, choose all the default options and let the installation finish.

## Building ZPIC
### Open the Command Line
Most Windows users - unlike Unix - never use the 'command line'. The 'comamnd line' is just another way to control the computer that uses text/keyboard commands only - rather then using the usually point-click mouse Windows interface. To start using the 'command line, go to the start menu. Scroll to the the folder entry starting with 'Visual Studio'. A screenshot to help:

![](https://github.com/atableman/zpic/blob/windows_test/doc/windows/win1.png)

 Click on this. Then click on entry labeled _Developer Command Prompt For..._. Another screenshot to help:
![](https://github.com/atableman/zpic/blob/windows_test/doc/windows/win2.png)

A black window should appear that looks like:
![](https://github.com/atableman/zpic/blob/windows_test/doc/windows/win3.png)

### Go To ZPIC Folder
Now, using the the _cd_ command, navigate to the folder where you placed the zpic code. In the screenshot below, ZPIC was in the folder `\lab_code\zpic` and the command `cd \lab_code\zpic` was used to change into this folder. Now Choose a version of ZPIC to compile - for example _em1d_ and go into the appropriate folder (using `cd em1d`)

### Make bin folder
Now, we mill make a new folder that will hold our built version of ZPIC. Let's call it _bin_. Make this folder by typing `mkdir bin` then go into this new folder using `cd bin`.

### Run CMake
Next, we will run the _CMake_ tool. That's easy:
```
cmake .. -G "NMake Makefiles"
```

Here's an example screen shot showing the previous commands:
![](https://github.com/atableman/zpic/blob/windows_test/doc/windows/win4.png)

### Run Make
After this completes, simply type:
```
nmake
```
### Done!
ZPIC is now ready use use! To Run ZPI, simple type:
```
zpic
```
### Editing source code.
Note: The CMake step only has to be done ONCE. When you make changes to the source code and want to rebuild ZPIC, simply run
```
nmake
```
again to rebuild ZPIC.
