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

where ``VER`` is the current version of `MELD <https://github.com/maccallumlab/meld-pkg>`_, and ``CUDA_VER`` is the installed version of cuda which is currently one of ``'75'``, ``'80'``, ``'90'``, ``'92'``, or ``'100'``. This will install MELD and all of its dependencies into a dedicated conda environment. Instructions on installing MELD from scratch can be found at `github <https://github.com/maccallumlab/meld>`_. 

The ideology behind MELD
=========================

The main idea behind MELD is to incorporate uncertain, ambiguous, or noisy experimental data into a simulation as distance restraints using Bayesian logic which is rooted in Bayes theorm:
::
	p(x|D) = p(D|X)p(x) ~ p(D|x)p(x)

where the likelihood of the data given a structure ``(p(D|x))`` and the prior ``p(x)``, which expresses our belief in the probability of a structure given no additional data, allow us to infer the posterior ``(P(x|D))``, which is the ensemble of likely structures given some data D. In our case, the prior is the Boltzmann distribution produced by the Amber forcefield. The likelihood is more complicated. 
In order to account for incorrect distance restraints from experimental data (whether due to ambiguity, noise, or uncertainty), MELD selects a subset of the distance restraints (a user-specified percentage, `N`) that are to be considered correct. For each structure ``x`` we assume that the `N` lowest energy restraints are most likely to be correct and are the most well-satisfied. The rest are ignored. 

By correctly sampling from the posterior distribution, ``p(x|D)``, MELD simultaneously infers the globally most likely configurations and the corresponding most likely restraints. MELD generates the minimum-free-energy ensemble by selecting the correct interpretation from the sea of possibilities.
A more complete explanation can be found `here <https://www.pnas.org/content/112/22/6985>`_. 


A basic example
================
This example will go through a simple MELD script that creates a system for a protein in an implicit water box in the absence of any external restraints. All of the files can be found in our github.

The setup script
----------------
The first step of any MELD simulation is to create a system using a setup.py script. The header contains our necessary imports:
::
    import meld
    import simtk.openmm as mm

We then set the total number of steps and the size of blocks we want our data to be saved in as global variables, as well as the number of replicas:
::
    N_STEP = 50000
    BLOCK_SIZE = 100
    N_REPLICAS = 32

The bulk of the work is done in a ``setup_system()`` function. The first step is to create the system. This can be done either from sequence
::
    protein_sequence = parse.get_sequence_from_AA1(filename='protein.dat')
    protein = system.ProteinMoleculeFromSequence(protein_sequence)

or from a PDB file
::
    protein = system.ProteinMoleculeFromPDB('protein.pdb')

Once all files are loaded we actually build the system:
::
    b = system.SytstemBuilder()
    s = b.build_system_from_molecules(protein, leap_header_cmds="source leaprc.water.tip3p")

This creates a system object of our protein in a TIP3P water box. MELD makes use of scalers to scale the temperature and the restraint energy along the replica ladder:
::
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 1.0, 300., 550.)

In this example we have no external restraints, so the temperature is scaled from 300 to 550 K along the whole of the ladder. Next we add in the options:
::
    options = system.RunOptions()
    options.implicit_solvent_model = 'obc'
    options.cutoff = 1.8
    options.timesteps = 25000
    options.minimize_steps = 5000

Here we have set our implict solvent model to OBC (a full list of available solvent models can be found here), our cutoff to 1.8, and our number of timesteps to 25000 with 5000 steps of minimization. We can calculate the overall length of our simulation:
::
    timestep x options.timesteps x N_STEPS x 10^-9 = microseconds
    2 fs x 25000 steps per block x 50000 blocks x 10^-9 = 2.5 microseconds

A 2 fs timestep is default. This can be changed using the option ``options.use_big_timestep = True``. This will make other changes, a full explanation of big timesteps can be found here. Next we create a data store:
::
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

MELD manages the replica exchange ladder using a runner, which must be created:
::
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1)
    remd_runner = master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS,ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)
 
For the most part these two steps will always be done like this. In order to communicate between replicas, we need to create and store a communicator using MPI:
::
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)
 
We then create and save the initial states of each replica:
::
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

The completely initialized data store is saved and the function returns the system: 
::
    store.save_data_store()
    return s.n_atoms

The last thing to do is call the ``setup_sytem`` function:
::
    setup_system()

Running our system
------------------
The first step is to load our conda environment:
::
    source activate meld-{VER}

We then setup the system:
::
    python setup.py

If we look into our directory we'll notice that a ``Data`` directory has been created which contains a Blocks directory that will store the system data during the simulation, a Backup directory, an dat files for the system, replica exchange runner, and the options. 
As MELD makes use of GPU's, the primary system requirement is a GPU-enabled system. We will first run through how to run MELD on a local machine before running on a cluster:
::
    launch_remd_multiplex --debug

This launches MELD using a single GPU. The debug option stores additional data that can be helpful in case of error. While running on a local machine may be useful for detecting issues with the system or for small systems using few (8 or fewer) replicas, it will be incredibly slow otherwise. Below is a sample script for running MELD on a cluster using the SLURM queueing system:
::
    #!/bin/bash
    #SBATCH --gres=gpu:4
    #SBATCH --nodes=8
    #SBATCH --ntasks=32
    #SBATCH --time=00-50:00
    #SBATCH -o protein.out
    #SBATCH -e protein.err
    #SBATCH --job-name="protein"
    
    module load cuda
    source activate meld-0.4.14
     
    mpirun -np 32 launch_remd --debug

The most important lines are the last two. We need to activate our meld environment on the clusters nodes. This example uses `mpirun` but `srun` is also valid depending on the queueuing system setup. Here we launch 32 processes/replicas, one per node. The .out and .err files contain the same data that was printed to screen on our local machine. 

Re-starting a simulation
------------------------
When running MELD on a cluster, the system may time out. We can restart a system from the last completed block:
::
    prepare_restart --prepare-run

Understanding MELD output
-------------------------
The bulk of the output is found in the Data directory. In our working directory however we find the ``remd.log`` as well as the ``.err`` and ``.out`` files, which are useful for identifying issues with the simulation. The trajectory of the 0th replica is located in the ``trajectory.pdb`` file in the Data directory. We can also generate a dcd file:
::
    extract_trajectory --extract_traj_dcd protein.dcd

This and other useful scripts are located in the `scripts` directory, which can be seen here.
