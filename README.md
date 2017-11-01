<!-- README.md file --> 
# Description

The **tkLayout-lite** represents a derived SW version of tracker geometry modelling tool and tracker performance analysis 
tool **tkLayout**, which has been fully developed within the CMS collaboration at CERN. Its main aim is to study 
performance of a new silicon tracker, replacing the current one in a view of the Phase-2 upgrade (early 2020's). In order 
to utilize such broad functionality of tkLayout SW in other experiments an intensive effort has been made to prepare a 
light version of the original tkLayout tool. Its main feature is a full modularity, experiment independence and 
comprehensive documentation.

The key features of tkLayout-lite may be described as follows: a tracker geometry is generated from a set of simple 
configuration files (cf. `GeometryManager` & `MainConfigHandler`), defining the module types, tracker layout and material 
composition. Various support structures are then automatically added and services routed to build a truly realistic 
tracking system (c.f. `Materialway`). The tkLayout-lite has been programmed with flexibility in mind to provide natural 
support for various detector geometries and a wide spectrum of analyses to be run with given layout and study the best 
trade-off among many figures of merit, such as: tracking resolution, material budget, power dissipation, cost etc. The 
tkLayout-lite SW is easy to be used thanks to its modular structure. The individual analysis modules, so-called analyzers, 
are run sequentially by so-called `AnalysisManager`, first to analyze the data and then to visualize the obtained results 
in a compact & user-friendly way. The visualization is provided via html format. Several web pages are generated and 
included on the main web-site. They may contain various explanation texts, pictures, graphs, tables etc. (cf. 
`RootWSite`, `RootWPage`, etc.). Several important analyzers to be explicitely mentioned: 

 * `AnalyzerResolution`: perform track resolution studies, which mathematically make use of a parabolic approximation in a 
 global chi2 fit technique to estimate the track parameters cov. matrix (parabolic approx. simplifies a more complex 
 approach to circle fitting in XY; a line is being used to simplify the track fitting in s-Z)
 * `AnalyzerPatternReco`: perform track propagation & estimation of tracker pattern recognition capabilities - makes use of 
 error propagation technique applied to a parabolic approximation. To estimate the background level (background hits) for 
 given tracker & accelerator a Fluka charged particles fluence map is assumed to be added to the tkLayout-lite on the input 
 -> not included within the standard SW config files
 * `AnalyzerGeometry`: report details about geometry layout
 * `AnalyzerMatBudget`: report details about expected tracker material budget 
 * `AnalyzerOccupancy`: report occupancy & date rates being estimated based on built-in geometry & Fluka charged particles 
 fluence map (assumed to be provided on the input to tkLayout-lite -> not included within the standard SW config files).   

Developed at CERN by Giovanni Bianchi, Nicoletta De Maio, Stefano Martina and Stefano Mersi. **tkLayout-lite** version 
developed at CERN by Zbynek Drasal in collaboration with the tkLayout team `https://github.com/tkLayout`.

# Getting the code
The code is accessible from an official github repository: `https://github.com/tkLayout` or author's github repository: 
`https://github.com/drasal/tkLayout` (use clone command to make a local copy):

    git clone -b masterLite https://github.com/tkLayout.git
    cd tkLayout

To join the effort on code development or to try out the newest features implemented during the tkLayout-lite development, 
clone directly the development branch `devLite` from author's github repository:

    git clone -b devLite  https://github.com/drasal/tkLayout.git
    cd tkLayout

# Before the compilation/run
Generally, one needs several libraries to be linked with the tkLayout-lite: a working version of **ROOT 6** (`root` and 
`root-config` should be in a user's path defined by `ROOTSYS` variable) and 2 specific **BOOST libraries**:

 * ROOT 6 library set (follow instructions on `https://root.cern.ch/downloading-root`)
 * `boost_filesystem`
 * `boost_regex`

On CERN lxplus machine, the procedure is quite straightforward, simply run a bash shell and source the following 
configuration file before the compilation:

    source setup_slc6.sh

On local Linux machine, the environmental variables need to be set first in a local `setup.sh` file:

 * `BOOST_LIB`: BOOST libs directory
 * `BOOST_INCLUDE`: BOOST include files directory
 * `BOOST_SUFFIX`: Usually no suffix on local machines, mainly used by complex computer systems to identify different versions of BOOST installations 
 * `ROOTSYS`: ROOT 6 directory (`root-config` used by CMake to configure ROOT libs and includes for compilation is assumed to be located in `ROOTSYS/bin`)

and then the file has to be sourced:

    source setup.sh

# Compilation/Install
Compilation using **CMake** (don't forget to setup variables using `setup.sh`/`setup_sl5.sh` beforehand):

    mkdir build     (all object files, help files, libs, ... will be kept here)
    cd build
    cmake ..        (generate makefile)
    make install    (or write make all + make install)
    make doc        (generate Doxygen-based documentation in doc directory)

    make uninstall  (if cleaning needed)
    rm *            (clean all content in build directory & restart if needed)

Create first a **build** directory in the tkLayout home directory. All make related content will be then created here. 
Change the working directory to **build** and call `cmake ..` command (directing CMake to the tkLayout uppermost directory 
so that it reads correctly the `CMakeLists.txt`, the main CMake configuration file). After all MakeFiles are generated by 
CMake, call `make all` (`make install`) to (re)compile the tkLayout-lite. Optional **install** will copy the executables 
into `tkLayout/bin` directory and create symbolic links in `${HOME}/bin` directory. Hence, if a user adds this directory 
into the `${PATH}` variable, the `tklayout` may be then executed from any directory on the system. To revert the 
**install** procedure, simply apply the **uninstall** option. Finally, in order to recompile the SW, run 
`make all`/`make install` command again.

## First-time install
If this is the first time that you install the tkLayout-lite, a few questions will be asked and a tkLayout configuration
file, `.tkgeometryrc`, will be created in a `$HOME` directory:

1. One needs to provide the destination directory for the html output of the program, such directory needs to be writable 
for the user and readable for the web server.
2. The set of momentum values to be used for the simulation
3. Project name
4. Author of generated results

## Documentation
Doxygen documentation is generated using **doc** option and all the generated html files are saved in `tkLayout/doc/html` 
directory.

# Geometry description
In order to ease the start-up with the **tkLayout-lite** SW, several reference geometry descriptions have been included 
within the SW package. They can be found in tkLayout-lite `geometries` directory, namely `geometries/CMS` (for CMS Ph2 
studies) and `geometries/FCC` (for FCC-hh studies). Create a `run` directory first, choose one of the geometries, e.g. 
`FCChh_v3.03`, copy its content to the `run` directory and run tkLayout-lite here (using the main geometry configuration 
file as an input parameter):

    cd tkLayout
    mkdir run
    cd run
    mkdir geometries
    cd ..
    cp -r geometries/FCC/FCChh_v3.03 run/geometries 
 
The geometry description is defined and structured across several configuration files, one of which is the main config 
file, e.g. `FCChh_v3.03.cfg`. The other files (like material description, individual modules description etc.) are included 
to the main file using the directive commands `@include`. In addition to the main geometry configuration file, there 
exist also a simulation run parameters file: `SimParms` and file specifying general settings: `.tkgeometryrc`. The 
first one is included to the main config file by directive `@include`, the latter is automatically created once 
tkLayout-lite is run for the first time in a user's home directory (see **First-time install** section). The 
`.tkgeometryrc` file may be deleted, tkLayout will automatically create a new one and ask the user to redefine all of its 
values. This file may also be directly modified by the user, the following parameters are particularly important:

    TKG_LAYOUTDIRECTORY="$userhome/tkLayout/run" -> Look for geometry layouts in this directory   
    TKG_STANDARDDIRECTORY="$userhome/tkLayout/run" -> Look for standard config files in this directory (xml, config, etc.)
    TKG_MOMENTA="1.0,2.0,5.0" -> Define individual values of transverse momenta in GeV/c
    TKG_PROJECT="FCC-hh" -> Define a project name

# Run tkLayout-lite
To learn about the tkLayout-lite functionality, simply run `tklayout` with no input parameters defined. The typical 
configuration being used to study a tracker geometry, material budget and tracker resolution is the following:

    cd run
    tklayout -n 1000 -N 1000 geometries/FCChh_v3.03/FCChh_v3.03.cfg 
  
The parameter `-n` defines statistics to be simulated for geometry calculations, the parameter `-N` defines statistics for 
material or resolution calculations.  

The simulation results are automatically saved in a `results` directory, which is either newly created or rewritten under 
geometry directory, e.g. `geometries/FCChh_v3.03/results`. Use any web browser on `results/index.html` file to display 
all results in a compact way. The content of such directory may be deleted as it's being automatically rewritten each time 
the tkLayout-lite SW is being run on the same geometry configuration. The way, how the web browser displays the final 
results depends on the cascade styles (being already predefined for the tkLayout SW). Hence, one needs to make a symbolic 
link to a **style** directory (predefined in tkLayout-lite main directory). Similarly, some general configurations 
(material database etc.) are predefined in **config** and **xml** directories. Make a symbolic link in a `geometries` 
directory to them.

# Summary of tkLayout-lite run options
The following tkLayout-lite functionality is currently supported:

    Analysis options:
    -h [ --help ]                       Display help info.
    --opt-file arg                      Specify an option file to parse the 
                                        program options from (in addition to 
                                        command line).
    -n [ --geometry-tracks ] arg (=100) Number of tracks for geometry 
                                        calculations.
    -N [ --material-tracks ] arg (=100) Number of tracks for material & 
                                        resolution calculations.
    -o [ --occupancy ]                  Report occupancy studies based on Fluka 
                                        data.
    -g [ --geometry ]                   Report geometry layout.
    -m [ --material ]                   Report material budget.
    -p [ --patternreco ]                Report pattern recognition capabilities 
                                        in given occupancy (using Fluka data).
    -r [ --resolution ]                 Report resolution studies.
    -a [ --all ]                        Report all studies.
    -e [ --extraction ] arg (=CMS FCC)  Extract tkLayout geometry to an XML file 
                                        to be used in CMS/FCC SW frameworks. 
                                        Supported values: CMS or FCC.
    --verbosity arg (=1)                Verbosity level (Overridden by the option
                                        'quiet').
    --quiet                             No output produced (Equivalent to 
                                        verbosity 0, overrides the 'verbosity' 
                                        option).
    --performance                       Outputs the CPU time needed for each 
                                        computing step (Overrides the option 
                                        'quiet').
 
In order to study the `--patternreco` and `--occupancy` options, a Fluka charged particles fluence map has to be provided 
in `config` directory. As the `--all` option is meant to study all available analysis options together, the Fluka map is 
also required there on input.

# Update
To get the latest development version (usually still under development or being tested) one simply types:

    git fetch
    git chechout masterLite (or devLite)
    make
    make install

# Compatibility disclaimer
The developers of **tkLayout-lite** and **tkLayout** SWs tried their best to provide the same core functionality in both 
"flavours" of tkLayout, nevertheless not all the analysis options are the same or being implemented. The tkLayout 
functionality is more focused on CMS Phase 2 upgrade tracker needs, the tkLayout-lite on the other hand on FCC-hh 
tracker needs.   
