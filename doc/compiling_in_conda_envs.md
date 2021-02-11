# Installing dependencies with conda

A convenient way to get the dependencies to compile and run MONC locally
is to use a [conda](https://docs.conda.io/en/latest/miniconda.html)
environment.

A conda environment which matches your local system can be conveniently
created with

```bash
python .conda/create_env.py
```

This command parses `.conda/environment.yml.meta` and adds the correct packages
for your system to an environment file in `.conda/environment.yml` and creates
a conda environment called `monc` which can be used for compiling MONC.
