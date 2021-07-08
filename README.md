**This github-based fork of MONC is rapidly changing, with automated testing,
docs, etc being set up**

[![Build Status](https://travis-ci.org/leifdenby/monc.svg?branch=master)](https://travis-ci.org/leifdenby/monc)

MONC is a highly scalable Large Eddy Simulation (LES) model that has been
developed to simulate clouds and turbulent flows at high resolution (~ 10s of
metres) on large domains. The MONC project was initially funded through a Joint
Weather and Climate Research Program (JWCRP) enabling post via NERC research
grant [NE/L01338X/1](http://gtr.rcuk.ac.uk/project/F1023ACF-76FF-4AD4-8E22-BA4A74451AD0),
which facilitated a collaboration between Met Office scientists and
computational scientists at [EPCC](http://www.epcc.ed.ac.uk/) to undertake the
development of MONC.

The main aim of the MONC project was to develop a high resolution atmospheric
process research model for the community that is user friendly and scalable on
modern high performance computing (HPC) systems. The scientific basis for MONC
is the Met Office Large Eddy Model (LEM) and the development involved the
complete re-write of the Met Office LEM using modern software design with
a flexible plug 'n play component based architecture with a focus on high
performance computing (HPC) scaling and efficiency.  '

# Compiling and Running on Leeds ARC machines

To make compilation more simple on the ARC systems, a script called `monc_compile_arc.sh` has been written which loads the necessary modules and compiles MONC using fcm. A second script, `submonc.sge`, can be used to submit a MONC task using the command `qsub submonc.sge` and similarly loads the required modules before submitting the job.

### MPI libraries and compilers on ARC

Testing has shown MONC to be incompatible with intel compilers at present (at the very least with intel 19.04), and so the compilation script uses the gnu compiler (version 4.8.5 by default on arc4, but tested and working with gcc 8.3.0 and gcc 6.3.0). 

MONC has been known to sometimes have issues with openmpi, and so mvapich2 is recommended. Mvapich2 is incapable of running with multiple threads however, and testing has shown successful runs using all three MPI variants available on arc4 (openmpi, intelmpi and mvapich2) so while mvapich2 is the recommended variant, its use remains only a recommendation rather than a necessity.

# Contributing

> **Note**: Github have written their own guides on [what git
is](https://guides.github.com/activities/hello-world/) and [forking
projects](https://guides.github.com/activities/forking/). A quick-start
guide for working specifically with MONC follows below, see those for more
detail.

To work on the MONC source code you'll first need your own copy:

1. [Create your own fork](https://github.com/Leeds-MONC/monc/fork) on
   github.com of the main Leeds-MONC repository at
   https://github.com/Leeds-MONC/monc/fork

2. Clone your fork locally

    ```bash
    $> git clone https://github.com/{your-github-username}/monc
    ```

Then you can start work on your new feature/bug fix :rocket: and
contribute it back to the main MONC repository:

3. Create your own branch or commit to `master` locally (remember you can
   do as many commits as you want, they are free and fast)

    ```bash
    $> git checkout -b cool-new-feature-for-monc
    (change files...)
    $> git add .
    $> git commit -m 'Add cool feature to MONC'
    ```

    ... (repeat `git add` and `git commit` as necessary)

4. Once your happy with your feature push it to github, check out and
   resolve any issues caught by the automated tests

    ```bash
    $> git push origin cool-new-feature-for-monc
    ```

5. [Create a pull request on
   github.com](https://github.com/Leeds-MONC/monc/pull/new/master) to
   request that your change is included into `master` on the main fork at
   https://github.com/Leeds-MONC/monc

6. Resolve any issues/ask for review by the rest of the Leeds-MONC
   community on your pull request. Once ready, **merge**!

7. Celebrate :smiley: :tada: :confetti_ball:, and thank you for your hard
   work! :star:

To make sure your fork is up-to-date with master (and `master` on your
laptop) you'll need to add the main MONC repository as a "remote" that git
knows about (the convention is to call this `upstream`) and "pull" down
any changes:

8. Add the main repository as an upstream (if you don't already have it):
    
    ```bash
    $> git remote add upstream https://github.com/Leeds-MONC/monc
    ```

9. Pull in the changes to your local `master` branch (or any other branch where
   you want recent changes into)

    ```bash
    $> git pull upstream master
    ```

10. To update `master` on your fork on github you can now push the changes to
    there

    ```bash
    $> git push origin master
    ```


# Keeping in sync with the Met Office Science Repository (MOSRS)

MOSRS uses SVN which can be used from git using `git svn`. The way this is done
is to set up a `svn-remote` (pointing to a specific SVN path, for example
`trunk`) and then fetching this content. Each commit on the MOSRS SVN will
then have an corresponding git commit. Adding and accessing `trunk` (or
any other branch) from MOSRS can be done with the script in
[misc/add_mosrs_remote.sh](misc/add_mosrs_remote.sh). This script ensures
that your local `master` is in sync with `master` on the Leeds fork,
fetches the remote content from MOSRS, creates a branch from that content
and ensures the history syncs up.

    ```bash
    $> ./misc/add_mosrs_remote.sh <mosrs_svn_remote_path> <local_branch_name>
    ```

You will need a user account for MOSRS to be able to run this command.
Once done the entire history on MOSRS will then be available on the to git
(and can be merged into existing branches etc).
