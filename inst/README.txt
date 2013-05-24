## Installing the GNU Scientific Library (GSL) on Windows. 

iBMQ employs functionality in GSL to improve the speed 
of various matrix and vector operations. Although it is
possible to compile GSL from source, we do not
recommend it. Precompiled binaries are available from

http://gnuwin32.sourceforge.net/packages/gsl.htm

After installing GSL, proper environment variables must be
established so that iBMQ can link to GSL.
In Windows 7 Ultimate, environment variables can be found in

Control Panels  > System and Security
		> System 
		> Advanced system settings 
		> Environment Variables

Note this path may vary by Windows version. If you are installing
the iBMQ package as a precompiled binary (.zip) file, you need
only edit the System variable "Path". To the path variable,
append /path/to/bin, where /path/to/bin leads to the GSL
(.dll) files. Note the following example "Path" variable.

System Variable "Path":

c:\Rtools\bin;
c:\Rtools\gcc-4.6.3\bin;
C:\Windows\system32;
C:\Windows;
C:\Windows\System32\Wbem;
C:\Windows\System32\WindowsPowerShell\v1.0\;
C:\Program Files\MiKTeX 2.9\miktex\bin\;
C:\Program Files\R\R-3.0.1\bin\i386;
C:\Program Files\GnuWin32\bin

In a sample installation, the last line 

"C:\Program Files\GnuWin32\bin"

contains the GSL binaries.

# Compiling iBMQ from source in Windows
Compiling iBMQ from source is only recommended for experienced users.
In addition to setting the "Path" System variable, one must also
establish User variables "GSL_INC" and "GSL_LIB" leading to 
GSL header and source files, respectively. 

In a sample installation:

GSL_INC: C:/Program\ Files/GnuWin32/include
GSL_LIB: C:/Program\ Files/GnuWin32/lib

Note the use of forward slashes separating directories. The character
"\ " handles spaces in directory names. 

After establishing environment variables for GSL, next one must
install the "Rtools" toolchain. This is available from 

http://cran.r-project.org/bin/windows/Rtools/

A required step in the Rtools installation is editting the "Path"
system variable, as described in the instructions

http://cran.r-project.org/bin/windows/Rtools/Rtools.txt

The sample "Path" variable above is an example of a working
Rtools "Path" variable, but this is machine and version dependent.

From a command prompt, navigate to the directory containing the
iBMQ source code. Execute "R CMD INSTALL" on the source code. This
will install iBMQ from source.


## Installing the GNU Scientific Library (GSL) on Unix-alikes

Visit http://www.gnu.org/software/gsl/ for instructions on how to
install GSL on Unix-alike platforms. iBMQ should
automatically detect the GSL libraries during compliation.

