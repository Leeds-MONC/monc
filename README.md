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


# Keeping in sync with the Met Office Science Repository (MOSRS)

MOSRS uses SVN which can be used from git using `git svn`. The way this is done
is to set up a `svn-remote` (pointing to a specific SVN path, for example
`trunk`) and then fetching this content. Each commit on the MOSRS SVN will then
have an corresponding git commit. Adding and accessing `trunk` from MOSRS can
be done with the following commands:

```bash
export LOCAL_NAME="mosrs-trunk"
git config --add svn-remote.${LOCAL_NAME}.url https://code.metoffice.gov.uk/svn/monc/main/trunk
git config --add svn-remote.${LOCAL_NAME}.fetch :refs/remotes/mosrs/trunk
git svn fetch ${LOCAL_NAME}
git checkout remotes/mosrs/trunk -b ${LOCAL_NAME}
unset LOCAL_NAME
```

You will need a user account for MOSRS to be able to run this command. Once
done the entire history of `trunk` on MOSRS will then be available on the git
branch called `mosrs-trunk`, which can then be merged into the current branch.
