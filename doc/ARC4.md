#
## A field guide to MONC with CASIM and SOCRATES on ARC4

Chris Symonds, Mark Richardson, Steven Boeing, Craig Poku and Leif Denby

23/4/2021

Aim: instructions on how to retrieve, compile and run MONC with CASIM and SOCRATES, and notes on internal workings of the coupling between MONC and SOCRATES/CASIM

# _Overview_

**CASIM** provides a two-moment(?) microphysics parameterisation which is called on every timestep to predict the formation and removal of water condensate (ice and liquid).

**SOCRATES** provides functionality to compute column-wise longwave absorption and emission of radiation represented by cooling and heating at every grid point.

The guide is split into three steps: a) retrieving a copy of MONC, b) compiling MONC on ARC4 and c) compiling MONC with SOCRATES and CASIM.

# a. Retrieving MONC

To get a copy of MONC from the Leeds fork on github is recommended you first [create your own fork](https://github.com/Leeds-MONC/monc/fork) and then clone that fork locally onto the computer where you are working:

```
$> git clone https://github.com/<your-github-username>/monc/
```


# b. Compiling and running MONC on ARC4

To compile and run MONC on ARC4 you will need to ensure that the correct versions of required libraries are loaded and then compilation to take place. This should all be taken care of by the script in [utils/arc/compile_monc.sh](../utils/arc/compile_monc.sh) which can be run with

```bash
$> bash utils/arc/compile_monc.sh
```

from the root of the repository. For completeness the steps the contents of that script are detailed below.

## 1. Make sure the dependencies of MONC are available

Issue the following commands to load the correct modules with `module`:

```bash
$> module purge
$> module load user
$> module switch intel gnu/4.4.7
$> module switch openmpi mvapich2
$> module load netcdf hdf5 fftw fcm
```

At time of writing CP compiled with the following versions: gnu compiler `4.4.7`, mvapich2 `2.3.1`, fcm `2019.09.0`, hdf5 `1.8.21`, netcdf `4.6.3` and fftw `3.3.8`

Notes on versions:

- gnu version: At time of writing Craig Poku's copy of MONC only works with `gnu/4.4.7`, but changes for `gnu/8.3.0`
- MPI implementation: `openmpi` or mvapich2 doesn't allow for multithreading, but Rachel Stratton has found MONC to be more stable with `mvapich` on ARC4
- fftw: due to licensing issues use of fftw is no longer the default Fourier transform used on the head of MOSRS `trunk`, instead ffte is used in current head-of-trunk on MOSRS (ffte is included with the MONC sourcecode on MOSRS `trunk`). In principle there is no issue of using FFTW and if your research can work with the licensing and it is installed you might be happier using it. MOSRS was not clearly labelled before they made that change.

## 2. Compile your copy of MONC

`fcm` is here wrapping `make`, while including correct .cfg-files in `fcm-make/` define set up the compilation environment

```bash
$> fcm make -f fcm-make/monc-arc4-gnu.cfg -N --ignore-lock –j4
```

## 3. Submitting a MONC run

Job submission script (and example is given in [../utils/arc4/submonc.sge](submonc.sge)) should have the following important features:

1. Job _walltime_ longer than _walltime_ in monc configuration
2. Module loads (or later `module load monc-dependencies`)
3. MVAPICH variables:
```
MONC_THREAD_MULTIPLE=0 # to do with "thread-pooling" in MONC
MV2_ENABLE_AFFINITY=0
MV2_SHOW_CPU_BINDING=1
MV2_USE_THREAD_WARNING=0
export MONC_THREAD_MULTIPLE MV2_ENABLE_AFFINITY \
MV2_SHOW_CPU_BINDING MV2_USE_THREAD_WARNING
```

The you execute your job submission script
```bash
$> qsub <you-job-script.sub>
```

NOTE: before submitting the job ensure that the "standard out" path is cleared, so that the job isn't restarted from a previous run.

## 4. Check that job is queued/running

You can now check that your job is running (the the "cluster" _group_ on ARC4 is called `c`):

```bash
$> qstat -g c
```

If this run completed succesfully you can now continue onto compiling and running MONC with SOCRATES and CASIM.

# c. Compiling MONC with SOCRATES and CASIM

Due to the license of SOCRATES and CASIM the sourcecode for both resides on the MOSRS (Met Office Science Repository) for which you will need access to using MONC with CASIM/SOCRATES. Once you have your credentials you can follow the steps below to set up `fcm` so that you can retrieve SOCRATES/CASIM with `fcm` and compile MONC with either or both componeents.

## 1. Setup and check SVN connection
Add MOSRS to `svn` (subversion) list of servers by adding the following lines to `~/.subversion/servers`

```
[groups]
metofficesharedrepos = code*.metoffice.gov.uk

[metofficesharedrepos]
username = <your-username>
store-plaintext-passwords=no

```

[https://code.metoffice.gov.uk/trac/monc/wiki/MoncDoc/MoncUserguide/MosrsSetup](https://code.metoffice.gov.uk/trac/monc/wiki/MoncDoc/MoncUserguide/MosrsSetup)

Now check the SVN connection (may be asked to cache the password, Craig found it best not to cache it encrypted)

```bash
$> svn info https://code.metoffice.gov.uk/svn/test
```

Caching password is covered in above MOSRS link and here: [http://cms.ncas.ac.uk/wiki/MonsoonSshAgent](http://cms.ncas.ac.uk/wiki/MonsoonSshAgent).

## 2. Check fcm keywords

[https://code.metoffice.gov.uk/trac/monc/wiki/MoncDoc/MoncUserguide/FcmKeyWords](https://code.metoffice.gov.uk/trac/monc/wiki/MoncDoc/MoncUserguide/FcmKeyWords)

Check what keywords are there by doing following command:

```bash
$> fcm keyword-print
```

If `socrates` or `casim` or not listed in these keywords you will to point `fcm` to the correct set of keywords. add the following in a file at `~/.metomi/fcm/keyword.cfg`:

```
location{primary, type:svn}[monc.x] = https://code.metoffice.gov.uk/svn/monc/main
browser.loc-tmpl[monc.x] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[monc.x] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /\*(.\*) \z)

location{primary}[casim.x] = https://code.metoffice.gov.uk/svn/monc/casim 
location{primary}[monc-doc.x] = https://code.metoffice.gov.uk/svn/monc/doc 
location{primary}[monc-postproc.x] = https://code.metoffice.gov.uk/svn/monc/postproc 
location{primary}[monc-scripts.x] = https://code.metoffice.gov.uk/svn/monc/scripts 
location{primary}[monc.x] = https://code.metoffice.gov.uk/svn/monc/main 
location{primary}[socrates.x] = https://code.metoffice.gov.uk/svn/socrates/main 
location{primary}[socrates.xm] = file:///home/d04/fcm/srv/svn/socrates.xm/main
``` 

The run `fcm keyword-print` again to check that the keywords now include `casim` and `socrates`.

## 3. Obtaining source-code for SOCRATES/CASIM and compiling MONC with these

In addition to the steps above, a copy of the source code for SOCRATES/CASIM is needed to compile MONC with CASIM/SOCRATES (both are treated identically in the comments that follow). The location of CASIM/SOCRATES is given through fcm configuration file (.cfg-file) which if provided through the command line (with –f) to fcm. This location can either be: a) a local filesystem path or b) a remote SVN path on MOSRS (with optional revision number from that SVN path to use). If option b) is used fcm will both check out the SOCRATES/CASIM code from MOSRS and compile with MONC

When calling `fcm` you need to include the `casim.cfg`, `socrates.cfg` or `casim_socrates.cfg` .cfg-files (if you're compiling with either CASIM, SOCRATES or both). Here we're including both CASIM and SOCRATES:

```bash
$> fcm make -f fcm-make/monc-arc4-gnu.cfg -f fcm-make/casim_socrates.cfg -N --ignore-lock –j4
```

You need to place the `casim`/`socrates`/`casim_socrates` fcm make config file _after_ the MONC file. If building with both casim and socrates, it needs a single fcm make configuration that deals with both rather than using two different config files.

The fcm configuration file does a number of things:

1. Instruct fcm to _extract_ code for CASIM and SOCRATES (makes sure the source files arencluded)
2. Include interaces for CASIM and SOCRATES, instead of place-holder routines
3. Define locations of where SOCRATES/CASIM is coming from
4. Set environment variables needed for CASIM/SOCRATES when being compiled to run inside MONC


## 2. Running MONC with CASIM and SOCRATES

To run MONC with CASIM and SOCRATES three things are needed in the model configuration (`.mcf`) file:

1. Flags to enable CASIM/SOCRATES and disable the functionality they replace

	```
	# required flags for CASIM
	simplecloud_enabled=.false.
	casim_enabled=.false.
	# required flags for SOCRATES:
	socrates_couple_enabled=.true.
	lwrad_exponential_enabled=.false. # turn off "bulk-calculation of radiation"
	```

2. Test-case specific parameters. These will define what microphysics processes to include specific configuration parameters for CASIM.

	```
	# number of scalar fields allocated, this needs to match the
 	# total number of tracers _required_. water vapour, graupel, etc
	number_q_fields=9
	```

3. External files providing reference profiles for radiation calculations (SOCRATES) and microphysics (CASIM). Craig Poku noted that providing _relative paths_ for these files based on where CASIM/SOCRATES is checked out within the MONC source tree (when fcm is used to fetch CASIM/SOCREATES sources from MOSRS) doesn't work and these configuration file parameters should point to where one has checked out the SOCRATES/CASIM source by hand.

## 3. Looking at the output

Can be done through ncview or similar. The outputs are in folder `diagnostic_files` and defined in the mcf file.