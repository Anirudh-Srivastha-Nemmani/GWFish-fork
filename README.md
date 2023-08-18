# GWFish

## About

This branch contains the edits done for using TEOBResumS in GWFish repository.

Author of this branch - Anirudh Srivastha Nemmani

## Instructions to install TEOBResumS environment

Detailed guide to create an environment to generate and run eccentric module of TEOBResumS.

- First create a conda environment using the following command conda create -n environment_name gsl , replace environment_name with whatever name you want and type y and press enter when asked. Note - Don't change the gsl part and make sure you have space between the environment_name and gsl .
- Now once the environment is created, activate it using conda activate environemnt_name .
- Now that you are in the environment, clone the gcc repository using the command git clone git://gcc.gnu.org/git/gcc.git
- Now go to the gcc directory using gcc . Now you have to checkout to the gcc-8 branch using git checkout origin/releases/gcc-8.
  - Note - gcc-8 is found to be working well with TEOBResumS, so this guide follows gcc-8, if you feel like you can make it work with higher versions of GCC, please feel free to update me.
- Before installing gcc, we need to check the pre-requisites using the command ./contrib/download_prerequisites
- Now go out of gcc directory and create a directory where you will build gcc.
- From the the build directory you have created in Step 6, run the following command ../gcc/configure --prefix=env_path --enable-languages=c,c++,fortran,go --disable-multilib . Note - env_path is the path to your environment. It would be in your anaconda envs folder. It would looking something similar to this /home/[your pc/account name]/anaconda3/envs/environment_name
- Once it's completed run make -j ncpus. Here ncpus is the number of cpus you want to use to run the make file. This will take a considerable amount of time, so it's recommended to use 4-8 cores. (Note - Make sure there are no lines that start with ***, *** indicates that the line is showing error output).
- After it's completed, run the command make install.
- Now run make clean to clean the build files.
- Now try to run gcc --version to check the installed version.
- Now exit the directory and install the libc libraries git repository using git clone git@github.com:hyperrealm/libconfig.git
- Now enter the git directory using cd libconfig/
- To produce configure file use the command autoreconf
- Now run the configure file using the command ./configure --prefix=env_path . Note - env_path is the path to your environment. It would be in your anaconda envs folder. It would looking something similar to this /home/[your pc/account name]/anaconda3/envs/environment_name
- Repeat step 8, 9, 10 in this folder.
- Now add this export in your bashrc file  export LD_LIBRARY_PATH=[path_to_lib_in_environment]:$PATH (you can access your bashrc using the following command vi ~/.bashrc , you can also use gedit instead of vim). [path_to_lib_in_environment] would to the path to the library in your anaconda envs folder. example /home/[your pc/account name]/anaconda3/envs/environment_name/lib
- Now save the file and source your bashrc using source ~/.bashrc.
- Now you are done installing all the requirements for teobresums eccentric branch.
- Now clone the teobresums repository using the command git clone https://bitbucket.org/eob_ihes/teobresums.git and enter the directory using cd teobresums
- Now checkout to the eccentric branch using the command git checkout eccentric
- To install the package it is best to refer to the instructions given in the teobresums installation guide as it would be updated from time to time. It would be summarised in the following steps.
- Enter into the python folder using the command cd Python. Now to install the cython module, run the following command python TEOBResumSWrap_setup.py build_ext --inplace . Make sure you have the .so file generated.
- Now open python shell using the command python3 and run the command import EOBRun_module , if it's imported properly then the eccentric package has been installed properly. You can also check it by running the Test.py in the examples folder. Note - when you are importing EOBRun_module you need to give the path of the directory where the cython module .so is present using the command sys.path.append('path_to_cython_module').
- If there are any issues in between the steps, contact me or Apratim.


<p align="center">
  <img src="gwfish-1.png" width="200" title="Logo">
</p>
Simulation of gravitational-wave detector networks with Fisher-matrix PE

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

Read the documentation [here](https://gwfish.readthedocs.io)!

## Citation

Please cite [GWFish publication](https://doi.org/10.1016/j.ascom.2022.100671) if you make use of the code:
```
@ARTICLE{DupletsaHarms2023,
        author = {{Dupletsa}, U. and {Harms}, J. and {Banerjee}, B. and {Branchesi}, M. and {Goncharov}, B. and {Maselli}, A. and {Oliveira}, A.~C.~S. and {Ronchini}, S. and {Tissino}, J.},
        title = "{GWFISH: A simulation software to evaluate parameter-estimation capabilities of gravitational-wave detector networks}",
        journal = {Astronomy and Computing},
        keywords = {General Relativity and Quantum Cosmology},
        year = 2023,
        month = jan,
        volume = {42},
        eid = {100671},
        pages = {100671},
        doi = {10.1016/j.ascom.2022.100671},
        archivePrefix = {arXiv},
        eprint = {2205.02499},
        primaryClass = {gr-qc},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2023A&C....4200671D},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Open GWFish tutorial in Google Colab

The tutorial notebook can be opend in Google Colab without the need to download locally any package. Here is the link: [notebook GWFish](<https://colab.research.google.com/github/janosch314/GWFish/blob/main/gwfish_tutorial.ipynb>)


This work is licensed under a [Creative Commons Attribution 4.0 International
License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
