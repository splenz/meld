==========================
Getting started with MELD
==========================

Installing MELD
================

The preferred way to install MELD is using `anaconda <https://www.anaconda.com/distribution/>`_ for python 3:
::
	conda config --add channels maccallum_lab omnia
	conda create -n meld-{VER}
	conda activate meld-{VER}
	conda install meld-cuda{CUDA_VER}

where ``VER`` is the current version of `MELD<https://github.com/maccallumlab/meld-pkg>`_, and ``CUDA_VER`` is the installed version of cuda which is currently one of ``'75'``, ``'80'``, ``'90'``, ``'92'``, or ``'100'``. This will install MELD and all of its dependencies into a dedicated conda environment. Instructions on installing MELD from scratch can be found at `github <https://github.com/maccallumlab/meld>`_. 

The ideology behind MELD
========================

The main idea behind MELD is to incorporate uncertain, ambiguous, or noisy experimental data into a simulation as distance restraints using Bayesian logic which is rooted in Bayes theorm .. math::
	p(x|D) = p(D|X)p(x) ~ p(D|x)p(x)

where the likelihood of the data given a structure ``(p(D|x))`` and the prior ``p(x)``, which expresses our belief in the probability of a structure given no additional data, allow us to infer the posterior ``(P(x|D))``, which is the ensemble of likely structures given some data D. In our case, the prior is the Boltzmann distribution produced by the Amber forcefield. The likelihood is more complicated. 
In order to account for incorrect distance restraints from experimental data (whether due to ambiguity, noise, or uncertainty), MELD selects a subset of the distance restraints (a user-specified percentage, `N`) that are to be considered correct. For each structure ``x`` we assume that the `N` lowest energy restraints are most likely to be correct and are the most well-satisfied. The rest are ignored. 

By correctly sampling from the posterior distribution, ``p(x|D)``, MELD simultaneously infers the globally most likely configurations and the corresponding most likely restraints. MELD generates the minimum-free-energy ensemble by selecting the correct interpretation from the sea of possibilities.
A more complete explanation can be found `here <https://www.pnas.org/content/112/22/6985>`_. 


A basic example
================
The first step of any MELD simulation is to create a system using a setup.py script.
