
Pyger Documentation
===================


This documentation is implemented as an ipython notebook. Please download and run this notebook and explore this functionality to get a feel for how Pyger works. Please visit http://ipython.org/notebook.html for more information on Python notebooks.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Pyger is a tool for calculating allowed hot dwell times and associated
cooldown times given spacecraft thermal constraints. More specifically
Pyger:

-  Uses a Monte-Carlo approach to sample a realistic ensemble of
   simulation starting temperatures.

-  Generates both typically achievable (50%) and best-case (90%) hot
   dwell times.

-  Generates typically achievable cooldown times for individual MSIDs.

-  Includes all xija-based thermal models and the ``pline`` mode for
   available hot time (not cooldown time).

-  Hot constraints implemented as Python classes derived from a base
   ConstraintModel class.

-  Cooldown constraints implemented as Python Numpy Arrays

-  Outputs data to pickle files for later use.

Pyger can be used in two complementary ways:

-  Command line usage: generate a single constraint plot for given set
   of conditions

-  Python scripting: automate the generation of many plots

Command line usage will be covered first, since the first step to create
the simulation inputs should to be executed at the command line. Then
the scripting interface will be discussed along with a the underlying
methods used to generate both maximum hot dwell and associated cooldown
times.

Command line usage
------------------


Some basic functionality can be performed directly from the command
line. This is useful for generating the inputs required for all Pyger
calculations, as well as testing to ensure Pyger is functioning
correctly. At this time, more detailed analyses are best performed
either within the ipython environment or via Python scripting.

The command line usage section cannot be executed directly from this
ipython notebook. Please execute these commands using a terminal in your
working directory. You will at least need to execute the ``pyger make``
command to generate the required simulation inputs, as described below,
which will be needed later for the pyger scripting section.

Pyger is available from within the Ska runtime environment. To
initialize the environment enter:

``% source /proj/sot/ska/bin/ska_envs.csh``

After this confirm that the pyger tool is available by asking for help:

::

    % pyger --help
    usage: pyger [-h] {make,sim} ...

    optional arguments:
      -h, --help  show this help message and exit

    Required sub-command (use exactly one):
      {make,sim}  Use "pyger <sub_command> --help" for help on sub-commands
        make      Make simulation inputs
        sim       Simulate constraints

This help message illustrates that the pyger command line tool has two
distinct functionalities:

-  Make the simulation inputs that are required for subsequent
   constraint simulations

-  Run a constraint simulation



Make simulation inputs
~~~~~~~~~~~~~~~~~~~~~~


In order to generate realistic set of dwell times pyger requires a
realistic set of starting conditions that are appropriate to the time of
the simulation. Since there is no longer a single "best-case" starting
temperature as there was for a 1-DOF TEPHIN model, pyger does this in
the following fashion:

-  Starting condition is a vector of temperatures (e.g. TEPHIN,
   TCYLAFT6, ...) or accumulated warm time

-  Need an ensemble of starting conditions based on realistic observing
   profiles leading up to start of simulation

-  Use a random sampling of actual profiles from the last year to
   generate propagated conditions:

-  Choose actual observation dwells to set conditions just prior to
   simulation start (maneuvers are filtered out of initial sampling set)

-  A default minimum dwell time of 1ksec is used to generate sampling
   set of final dwells for propagation profiles

-  For each profile start propagation 3 days before simulation start
   time

-  Start with as-observed temperatures at the beginning of the 3 day
   profile

-  Use thermal model to propagate temperatures during 3 day profile at
   simulation start date/time so simulation starts with representative
   conditions

-  This works for any thermal model including PSMC with variable
   CCD-count dependence

Actually generating the simulation inputs file is done with the command
``pyger make``. First look at the available options:

::

    % pyger make --help
    usage: pyger make [-h] [--sim-file SIM_FILE] [--start START] [--stop STOP] [--n-days N_DAYS]

    optional arguments:
      -h, --help           
                Show this help message and exit
      --sim-file SIM_FILE      
                Output file (default = sim_inputs.pkl)
      --start START       
                Sim input start time (default = Now - 1 year)
      --stop STOP         
                Sim input stop time (default = Now - 10 days)
      --n-days N_DAYS      
                Number of days to propagate prior to perigee exit (default = 3)
      --min-dwell-sec      
                Minimum number of seconds to use as a propagation ending dwell (default = 1000)
      --max-dwell-num      
                Maximum number of sample dwells to save (default = 300)

Now run it, which takes about two minutes using the recommended default
settings:

::

    % pyger make
    Finding suitable propagation ending dwells between:
      Start:2012:307:18:25:17.309
      Stop:2013:306:18:25:17.309
    Limiting number of sample dwells to 300
    Found 300 suitable propagation ending dwells at least 1000 seconds long.
    Assembling simulation propagation data for model: pline
      Fetching state values: pitch
      Fetching telemetry for: aosares1
    Assembling simulation propagation data for model: tank
      Fetching state values: pitch
      Fetching telemetry for: pftank2t
    Assembling simulation propagation data for model: minus_yz
      Fetching state values: pitch
      Fetching telemetry for: pmtank3t, tmzp_my, tephin, tcylaft6
    Assembling simulation propagation data for model: psmc
      Fetching state values: pitch, ccd_count, fep_count, clocking, vid_board, simpos
      Fetching telemetry for: 1pdeaat
    Assembling simulation propagation data for model: dpa
      Fetching state values: pitch, ccd_count, fep_count, clocking, vid_board, simpos
      Fetching telemetry for: 1dpamzt
    Assembling simulation propagation data for model: minus_z
      Fetching state values: pitch
      Fetching telemetry for: tcylaft6, tcylfmzm, tephin, tfssbkt1, tmzp_my
    Writing simulation inputs to sim_inputs.pkl

Note that the ``--start`` and ``--stop`` parameters can be specified in
any valid ``DateTime`` format. Occasionally you will receive an out of
bounds error at this step. If you receive this error, just reissue the
``pyger make`` command. This is a related to the random manner in which
dwells are chosen.

Run Constraint Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~


From the command line, constraint simulation and plot generation is done
with ``pyger sim``. First look at the available options:

::

    %pyger sim --help
    usage: pyger sim [-h] [--start START] [--n-sim N_SIM] [--dt DT]
                     [--max-tephin MAX_TEPHIN] [--max-tcylaft6 MAX_TCYLAFT6]
                     [--max-1pdeaat MAX_1PDEAAT] [--max-1dpamzt MAX_1DPAMZT]
                     [--max-pftank2t MAX_PFTANK2T] [--n-ccd N_CCD]
                     [--max-aacccdpt MAX_AACCCDPT]
                     [--max-dwell-ksec MAX_DWELL_KSEC] [--sim-file SIM_FILE]
                     [--plot-file PLOT_FILE] [--plot-title PLOT_TITLE]

    optional arguments:
      -h, --help                       Show this help message and exit
      --start START                    Simulation start time
      --n-sim N_SIM                    Number of simulation points
      --dt DT                          Delta time for sims (sec)
      --max-tephin MAX_TEPHIN          TEPHIN planning limit (degF)
      --max-tcylaft6 MAX_TCYLAFT6      TCYLAFT6 planning limit (degF)
      --max-1pdeaat MAX_1PDEAAT        1PDEAAT planning limit (degC)
      --max-1dpamzt MAX_1DPAMZT        1DPAMZT planning limit (degC)
      --max-pftank2t MAX_PFTANK2T      PFTANK2T planning limit (degF)
      --n-ccd N_CCD                    Number of ACIS CCDs
      --max-aacccdpt MAX_AACCCDPT      ACA CCD planning limit (degC)
      --max-dwell-ksec MAX_DWELL_KSEC  Max allowed dwell time (ksec)
      --sim-file SIM_FILE              Simulation inputs pickle file
      --plot-file PLOT_FILE            Output plot file name
      --plot-title PLOT_TITLE          Title on output plot

Then run it, which takes several minutes for the default settings:

::

    %pyger sim
    MINUS_YZ: calculating start temps for 300 dwells
    MINUS_YZ: simulating 500 dwells
    PSMC: calculating start temps for 300 dwells
    PSMC: simulating 500 dwells
    PLINE: calculating warm time for 300 dwells
    PLINE: simulating 500 dwells
    DPA: calculating start temps for 300 dwells
    DPA: simulating 500 dwells
    TANK: calculating start temps for 300 dwells
    TANK: simulating 500 dwells
    ACA: calculating start temps for 300 dwells
    ACA: simulating 500 dwells
    Writing constraint plot file constraints.png

Note that the ``--start`` parameter can be specified in any valid
``DateTime`` format.

Python Scripting
----------------


Hot Dwell Simulations
~~~~~~~~~~~~~~~~~~~~~


The more powerful way to use Pyger is via Python scripting and plotting.
This allows for the generation of complex plots in a reproducible way.

As for the command line interface one requires an initial simulation
inputs file.

The key function which is used to perform a constraint simulation is
``pyger.calc_constraints``.

This function has a lot of parameters but all of them have sensible
defaults so normally only a few need to be set. Here is a sample script
which is available on the HEAD and GRETA networks as
/proj/sot/ska/share/pyger/examples/plot\_by\_month.py:
In[1]:
.. code:: python

    """Plot best-case constraints for 2014 in 2-month intervals"""
    import matplotlib.pyplot as plt
    
    import pyger
    
    # Set up a few values that might be varied for different cases
    max_tephin = 154.0
    max_tcylaft6 = 105.0
    max_1pdeaat = 52.5
    max_1dpamzt = 33.0
    max_pftank2t = 93.0
    n_ccd = 6
    max_dwell_ksec = 200
    n_sim = 500   # Use 5000 for final plots but 500 for testing
    
    # Set up the font sizes and initialize figure
    plt.rc("axes", labelsize=10, titlesize=12)
    plt.rc("legend", fontsize=10)
    plt.figure(1, figsize=(6,4.5))
    plt.clf()
    
    # Step through 2014 in 2-month intervals
    for month in range(0, 12, 2):
        start = '2014-%02d-01T00:00:00' % (month+1)
        print '*** Calculating constraints for %s ***' % start
        constraints = pyger.calc_constraints(start=start,
                                             max_tephin=max_tephin,
                                             max_tcylaft6=max_tcylaft6,
                                             max_1pdeaat=max_1pdeaat,
                                             max_1dpamzt=max_1dpamzt,
                                             max_pftank2t=max_pftank2t,
                                             n_ccd=n_ccd,
                                             n_sim=n_sim,
                                             max_dwell_ksec=max_dwell_ksec)
        dwell1_stats = constraints['all'].dwell1_stats
        plt.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000, label=start[:7])
    
    # Finish making the plot with labels and then save to a file
    plt.title('Best-case dwells for 2014 by month')
    plt.legend(loc='upper center')
    plt.xlabel('Pitch (deg)')
    plt.ylabel('Dwell (ksec)')
    plt.grid()
    plt.ylim(0, max_dwell_ksec * 1.05)
    plt.savefig('dwells_2014_month.png')



.. parsed-literal::

    *** Calculating constraints for 2014-01-01T00:00:00 ***
    MINUS_YZ: calculating start temps for 300 dwells
    MINUS_YZ: simulating 500 dwells



.. image:: index_files/index_12_1.png


The output of the ``pyger.calc_constraints`` function needs a little
explanation. As can be seen from the script this function returns a
single value which is a Python dict containing a key/value pair for each
of the constraint models. Currently there are six constraint models
named ``minus_yz``, ``psmc``, ``dpa``, ``tank``, ``aacccdpt`` and
``pline``. The easiest way to understand the output is to run the
``plot_by_month.py`` script as shown and then start examining the
``constraints`` variable.
In[2]:
.. code:: python

    constraints





.. parsed-literal::
    {'aca': <pyger.pyger.ConstraintAca at 0x112a10d0>,
     'all': <pyger.base.ConstraintModel at 0x10803410>,
     'dpa': <pyger.pyger.ConstraintDPA at 0x112a1f90>,
     'minus_yz': <pyger.pyger.ConstraintMinusYZ at 0x11263350>,
     'pline': <pyger.pyger.ConstraintPline at 0x112a1a50>,
     'psmc': <pyger.pyger.ConstraintPSMC at 0x112a1c50>,
     'tank': <pyger.pyger.ConstraintTank at 0x112a1290>}



So there are six elements in the dictionary, each of which is a
ConstraintModel object. In addition to the six actual constraints there
is a special one called "all" which represents the merging of these
constraints, i.e. taking the shortest allowed dwell time among all
constraints. A ConstraintModel object has a number of attributes that
can be discovered with the ``dir`` built in python function.
In[3]:
.. code:: python

    dir(constraints['dpa'])





.. parsed-literal::
    ['__class__',
     '__delattr__',
     '__dict__',
     '__doc__',
     '__format__',
     '__getattribute__',
     '__hash__',
     '__init__',
     '__module__',
     '__new__',
     '__reduce__',
     '__reduce_ex__',
     '__repr__',
     '__setattr__',
     '__sizeof__',
     '__str__',
     '__subclasshook__',
     '__weakref__',
     '_calc_model',
     '_get_init_comps',
     '_get_model_temps',
     '_get_states1',
     '_set_init_comps',
     'calc_dwell1_T0s',
     'calc_dwell1_stats',
     'calc_dwells1',
     'calc_dwells2',
     'dwell1_stats',
     'dwells1',
     'find_limit_crossings',
     'find_long_hot_pitch_dwells',
     'limits',
     'max_dwell_ksec',
     'model_spec',
     'msids',
     'n_ccd',
     'n_sim',
     'name',
     'sim_inputs',
     'start',
     'state_cols',
     'times']


In[4]:
.. code:: python

    dir(constraints['all'])





.. parsed-literal::
    ['__class__',
     '__delattr__',
     '__dict__',
     '__doc__',
     '__format__',
     '__getattribute__',
     '__hash__',
     '__init__',
     '__module__',
     '__new__',
     '__reduce__',
     '__reduce_ex__',
     '__repr__',
     '__setattr__',
     '__sizeof__',
     '__str__',
     '__subclasshook__',
     '__weakref__',
     '_calc_model',
     '_get_model_temps',
     '_set_init_comps',
     'calc_dwell1_T0s',
     'calc_dwell1_stats',
     'calc_dwells1',
     'calc_dwells2',
     'dwell1_stats',
     'dwells1',
     'find_limit_crossings',
     'find_long_hot_pitch_dwells',
     'limits',
     'max_dwell_ksec',
     'n_sim',
     'name',
     'sim_inputs',
     'start']



The most important are:

-  dwells1: dwell pitch, duration, and constraint name for each
   simulated dwell

-  dwell1\_stats: binned statistics for the dwells

These are both NumPy record arrays (i.e. a table). The available columns
for each indivdual model are:
In[5]:
.. code:: python

    print constraints['dpa'].dwells1.dtype.names



.. parsed-literal::

    ('pitch', 'T0s', 'times', 'Ts', 'constraint_name', 'duration')

In[6]:
.. code:: python

    print constraints['dpa'].dwell1_stats.dtype.names



.. parsed-literal::

    ('pitch_bin0', 'pitch_bin1', 'pitch', 'dur50', 'dur90')


The available columns in ``constraints['all'].dwells1`` do not include
the ``T0s``, ``times`` and ``Ts`` columns since these latter columns
contain data specific to each simulation and are not relevant to the
composite maximum dwell time.
In[7]:
.. code:: python

    print constraints['all'].dwells1.dtype.names



.. parsed-literal::

    ('pitch', 'duration', 'constraint_name')


The columns are accessed as follows, for example:

::

    dwells1 = constraints['all'].dwells1
    plt.plot(dwells1['pitch'], dwells1['duration'], '.')

The code used by the pyger module to make a standard plot is
illustrative of useful techniques for scripting and plotting. This
function will be used below to generate a composite plot of hot dwell
capability using a new set of limits. Note the filtering by constraint
name:
In[8]:
.. code:: python

    def plot_dwells1(constraint, plot_title=None, plot_file=None, figure=1):
        """Make a simple plot of the dwells and dwell statistics for the constraint.
    
        :param constraint: ConstraintModel object (e.g. constraints['all'])
        :param plot_title: plot title
        :param plot_file: output file for plot
        :param figure: matplotlib figure ID (default=1)
        """
        plt.rc("axes", labelsize=10, titlesize=12)
        plt.rc("xtick", labelsize=10)
        plt.rc("ytick", labelsize=10)
        plt.rc("font", size=10)
        plt.rc("legend", fontsize=10)
    
        dwells1 = constraint.dwells1
        dwell1_stats = constraint.dwell1_stats
        plt.figure(figure, figsize=(6,4))
        plt.clf()
        names = ('none', '1pdeaat', 'tcylaft6', 'tephin', 'pline',
                 '1dpamzt', 'pftank2t')
        colors = ('r', 'g', 'k', 'c', 'b', 'm', 'y')
        for name, color in zip(names, colors):
            ok = dwells1['constraint_name'] == name
            plt.plot(dwells1['pitch'][ok], dwells1['duration'][ok] / 1000., '.' + color,
                     markersize=3, label=name, mec=color)
        plt.plot(dwell1_stats['pitch'], dwell1_stats['dur50'] / 1000., '-r')
        plt.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000., '-m')
        plt.grid()
        plt.title(plot_title or '')
        plt.xlabel('Pitch (deg)')
        plt.ylabel('Dwell (ksec)')
        plt.legend(loc='upper center')
        plt.ylim(constraint.max_dwell_ksec * -0.05, constraint.max_dwell_ksec * 1.05)
        plt.subplots_adjust(bottom=0.12)
        if plot_file:
            logger.info('Writing constraint plot file {0}'.format(plot_file))
            plt.savefig(plot_file)

Recalculate the constraints object using a new set of limits and store
it in the object named ``hot_constraints``.
In[9]:
.. code:: python

    models = ('minus_yz', 'aca', 'tank', 'psmc', 'dpa', 'pline')
    
    hot_constraints = pyger.calc_constraints(start='2014:001',
                                                       n_sim=1000,
                                                       dt=1000.,
                                                       max_tephin=154.0,
                                                       max_tcylaft6=105.0,
                                                       max_1pdeaat=52.5,
                                                       max_1dpamzt=33,
                                                       max_pftank2t=93.0,
                                                       max_aacccdpt=-15.0,
                                                       n_ccd=5,
                                                       sim_file='sim_inputs.pkl',
                                                       max_dwell_ksec=400.,
                                                       min_pitch=45,
                                                       max_pitch=169,
                                                       bin_pitch=2,
                                                       constraint_models=models)



.. parsed-literal::

    MINUS_YZ: calculating start temps for 300 dwells
    MINUS_YZ: simulating 1000 dwells


Now plot the composite of all constraints using the built in
``plot_dwells1`` function. This plot shows the Monte-Carlo sampling of
different starting conditions at each pitch. The different dot colors
represent the limiting constraint. The lower red line shows the median
allowed dwell time within each 2 degree pitch bin. The upper magenta
line shows the 90% value in each bin and represents a reasonable "best
case" at that pitch.
In[10]:
.. code:: python

    pyger.plot_dwells1(hot_constraints['all'])




.. image:: index_files/index_28_0.png


Cooldown Dwell Simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~


As hot dwell capability in all pitch regions becomes more constrained,
we need to start considering the required cooling time to adequately
balance hot and cool time for all relevant constraints.

Cooldown time is calculated on an individual constraint basis so only
the MSID of interest should have a limiting constraint for each set of
hot and cool dwell simulations when cooldown time needs to be
calculated. Mixing multiple constraints in this context will result in
unintuitive and/or unusable results.

Since we need to start with a clean set of hot dwell simulations that
contain only one limiting constraint, first create the
``hot_constraints`` object using the same syntax used above. We'll focus
on TCYLAFT6 only for this example by setting the limits for all other
constraints well beyond an observable value, 999.0 should be more than
adequate. When focusing only on one constraint, one should also limit
the models run to the model used to calculate the constraint of interest
using the ``constraint_models`` keyword as shown below.

Since generating Pyger data can be time consuming, a method for saving
both the heatup and cooldown data has been included. This tutorial will
show how to save this data. If a large survey of dwell capability is
being generated, such as for the development of an LTS, it is
recommended that the user generate individual MSID data for all msids,
dates, and limits separately. This data can then be read in and combined
according to the users' needs.
In[107]:
.. code:: python

    import matplotlib.pyplot as plt
    
    from Chandra.Time import DateTime
    
    import pyger
    
    hot_constraints = pyger.calc_constraints(start='2014:001',
                         n_sim=1000,
                         dt=1000.,
                         max_tephin=999.0,
                         max_tcylaft6=105.0,
                         max_1pdeaat=999.0,
                         max_1dpamzt=999.0,
                         max_pftank2t=999.0,
                         max_aacccdpt=999.0,
                         n_ccd=4,
                         sim_file='sim_inputs.pkl',
                         max_dwell_ksec=400.,
                         min_pitch=45,
                         max_pitch=169,
                         bin_pitch=2,
                         constraint_models=('minus_yz',))



.. parsed-literal::

    MINUS_YZ: calculating start temps for 300 dwells
    MINUS_YZ: simulating 1000 dwells


The ``save_pyger_pickle`` method saves all pickleable attributes of the
constraints object.
In[108]:
.. code:: python

    pyger.save_pyger_pickle(hot_constraints, 'demo_dwell1.pkl')

Since we only included one model in the ``constraints_models`` keyword,
only the ``minus_yz`` data is returned. The data stored by the ``all``
key will only include the impacts of limiting TCYLAFT6.
In[109]:
.. code:: python

    hot_constraints





.. parsed-literal::
    {'all': <pyger.base.ConstraintModel at 0x252e8290>,
     'minus_yz': <pyger.pyger.ConstraintMinusYZ at 0x220f9810>}



Below is a plot of max hot dwell duration using the built in plotting
routine discussed and shown above.
In[110]:
.. code:: python

    pyger.plot_dwells1(hot_constraints['minus_yz'])




.. image:: index_files/index_37_0.png


Generating the cooldown data is accomplished similarly to how the heatup
data is generated. Note that three objects are returned and each of
these objects are Numpy arrays, not ``ConstraintModel`` class types.

-  ``cool_constraints`` contains all of the relevant data from the
   cooldown simulations.

-  ``coolstats`` includes relevant cooldown dwell statistical data.

-  ``hotstats`` includes the relevant hot dwell durations for comparison
   to the cooldown data.

Although the hot dwell data is used to generate the starting conditions
for the cooldown simulations, the maximum hot dwell durations initially
calculated from these data are not directly relevant to the cooldown
simulations. The original maximum hot durations are general statistical
values (90th percentile line for example) for all pitch regions and do
not include the actual ending conditions necessary to serve as initial
conditions for the cooldown simulations. This will be discussed an
illustrated in detail below.

First lets dissect each of the keywords used by the
``pyger.calc_constraints2`` function before running it.

-  ``start='2013:001'``: This is the date for which the cooldown time is
   calculated and will usually match the start time used for the
   cooresponding hot dwell calculations. It is conceivable that one
   would want to adjust this date to account for the time used during
   the hot dwell calculations, however the difference in seasonal
   impacts between the start and end of the hot dwell should be
   negligible for dwell durations that the xija models are calibrated
   for.

-  ``n_ccd=None``: If this is set to ``None``, pyger uses the number of
   ACIS CCDs specified during the hot dwell simulations, if relevant to
   the model.

-  ``max_dwell_ksec=600``: This is the length of time in kiloseconds
   each cooldown simulation is run. Since it often takes longer to cool
   down a constraint compared to the time to reach the temperature
   limit, this keyword should be set to a time in kiloseconds much
   longer than the time used for the hot dwell; 600 kiloseconds is often
   adequate for the cooldown simulations.

-  ``pitch_num=50``: This is the number of pitch values simulated in the
   'not hot' regions for each hot dwell simulation. This results in a
   total of NxM cooldown simulations where N is the number of hot
   simulations and M is the number specified in ``pitch_num``.

-  ``pitch_range=None``: This is the range of pitch values within which
   to simulate cooldown dwells. When this keyword is set to ``None``,
   Pyger uses all pitch values previously used in the hot dwells
   simulations which did not reach the specified temperature limit. If
   this keyword is set to a two element tuple or list, then Pyger
   selects a uniform distribution of pitch values within the two pitch
   values in this tuple or list. This is useful in cases where the model
   configuration changes between the hot and cool dwell simulations,
   particularly if the heatup dwell uses 6 ACIS CCDs and the cooldown
   dwell uses 4 ACIS CCDs, or visa versa.

-  ``hot_dwell_temp_ratio=1.0``: This keyword governs the starting
   conditions used for each cooldown simulation. This ratio specifies
   the temperature ratio between the HOT dwell starting temperature
   (when cooler) and the HOT dwell ending temperature (at the
   temperature limit). A ratio of 1.0, as in this case, specifies that
   the model temperatures at the hot dwell ending temperature (for the
   MSID/constraint of interest) should be used for the initial
   conditions for the cooldown simulation. A ratio of 0.8 specifies that
   the starting conditions for the cooldown simulation should be taken
   when the MSID of interest reached a delta temperature ratio of 0.8
   between the starting and ending HOT dwell simulation. This can be
   useful if the user is interested in how the heatup and cooldown times
   compare when operating away from the temperature limit. Note that the
   ``hot_dwell_temp_ratio`` should not be less than
   ``1.0 - T_cool_ratio``.

-  ``T_cool_ratio=0.9``: This keyword governs when a cooldown simulation
   reaches the 'cooled' state. This ratio is also relative to the
   temperatures spanned during the hot dwell simulation. A ratio of 0.9
   means that the MSID, or constraint, of interest has cooled when the
   temperature reaches 90% of the distance between the temperature limit
   and the original starting temperature in the hot dwell simulation. It
   is recommended that a value of 1.0 not be used since Pyger will then
   only show which pitch regions cool to the very lowest temperatures
   and will show some cooling regions as neutral. The goal of these
   cooldown simulations is to show which regions cool in general; the
   regions that provide the most/fastest cooling will invariably show
   shorter cooling times regardless of what ``T_cool_ratio`` is used.

-  ``constraint_models``\ =('minus\_yz',): This is currently input as a
   tuple and should include only the model of interest. Note that the
   comma after the first entry within the parentheses is necessary,
   otherwise Python will not interpret this as a tuple.

-  ``msids=('tcylaft6',)``: This specifies the MSID, or constraint, of
   interest. As with the ``constraint_models``, this should be a one
   element tuple.

It is understood that the use of tuples for the last two keywords adds
unnecessary complexity; this construct will likely be simplified in the
future.

NOTE: There is a bug in how the ``hot_temp_dwell_ratio`` is used for
values less than 1.0, for now only use ratios of 1.0.

Next, take a look at the output of this function.
In[111]:
.. code:: python

    cool_constraints, coolstats, hotstats = pyger.calc_constraints2(
                                                        hot_constraints,
                                                        start='2013:001',
                                                        n_ccd=None,
                                                        max_dwell_ksec=600.,
                                                        pitch_num=50,
                                                        pitch_range=None,
                                                        hot_dwell_temp_ratio=1.0,
                                                        T_cool_ratio=0.9,
                                                        constraint_models=('minus_yz',),
                                                        msids=('tcylaft6',))



.. parsed-literal::

    tcylaft6 - MINUS_YZ: simulating 50 cooldown dwells for each of the 116 hot dwells
    Using a hot dwell temperature ratio of 1.0 to determine cooldown dwellstarting conditions


The output echos some of the relevant input information, and shows the
progress as each hot dwell is used to seed a set of cooldown
simulations. Although, due to the random nature of how the hot dwells
are selected, the actual number hot hot dwells used will change with
each run, this example should find roughly a hundred or so hot dwells to
use.

Next you will see that the first output, assigned to
``cool_constraints`` is actually a one element list. This is also an
example of unnecessary complexity that will be removed in the future.
For now, we'll just assign the sole element to the same variable name.
In fact, each of the outputs are single element lists so the single
element will be assigned the variable name of its parent.
In[112]:
.. code:: python

    cool_constraints.keys()





.. parsed-literal::
    ['tcylaft6']


In[113]:
.. code:: python

    cool_constraints = cool_constraints['tcylaft6']
    coolstats = coolstats['tcylaft6']
    hotstats = hotstats['tcylaft6']

Save these arrays to a pickle file directly using the pickle module
In[114]:
.. code:: python

    import pickle
    saved_data = {'msid':'tcylaft6',
                  'cool_constraints':cool_constraints,
                  'hotstats':hotstats,
                  'coolstats':coolstats}
    pickle.dump(saved_data, open('demo_dwell2.pkl', 'w'), protocol=2)

Next lets review the columns/fields available in the
``cool_constraints`` array. Each row corresponds to one hot dwell used
to seed a set of cooldown simulations. The columns are the names shown
below and correspond to the data stored in each element. Since there are
many cooldown simulations run for each hot simulation, some of these
individual elements will contain arrays.

-  ``dwell1_ind``: This is the index into the original ``dwells1`` array
   in the hot dwells constriaints object used as the seed for this set
   of cooldown dwells.

-  ``dwell1_pitch``: This is the original pitch used during the hot
   dwell.

-  ``msid``: This is the MSID, or constraint, of interest.

-  ``dwell1_duration_delta``: This is the duration during the HOT dwell
   that is equivalent to the time between the 'cooldown' temperature (as
   governed by the ``T_cool_ratio``) and the starting conditions for
   each the cooldown simulations. This returns a single value for each
   row.

-  ``dwell1_duration``: This returns the full duration, from the initial
   temperature to the limit, for the original hot dwell. This returns a
   single value for each row.

-  ``T_dwell1_0``: This contains the starting conditions for the initial
   hot dwell simulation. This returns a list of values, where each
   element corresponds to the initial temperature for each node in the
   model.

-  ``T_dwell2_0``: This contains the starting conditions for the
   cooldown dwell simulation. This returns a list of values, where each
   element corresponds ot hte initial temperature for each node in the
   model.

-  ``dwell2_pitch_set``: This is the set of 'not hot' pitch values used
   for each of the cooldown simulations.

-  ``dwell2_times``: This is the set of cooldown durations, where each
   duration represents the time it took to go from the cooldown starting
   conditions to the 'cooled' temperature (governed by the
   ``T_cool_ratio``). This returns a list of times, with a length equal
   to the number of cooldown simulations run. Also each of these times
   corresponds, in order, to the list of pitch values in
   ``dwell2_pitch_set``. Not all cooldown simulations will reach a
   'cooled' state. Simulations in neutral pitch regions will finish the
   simulation when the maximum dwell time length ``max_dwell_ksec`` has
   been reached.

-  ``dwell2_cool_temp``: This is the 'cooled' temperature that each of
   the cooldown simulations need to reach in order to be considered
   cooled.

These fields/columns are shown below.
In[115]:
.. code:: python

    cool_constraints.dtype.names





.. parsed-literal::
    ('dwell1_ind',
     'dwell1_pitch',
     'msid',
     'dwell1_duration_delta',
     'dwell1_duration',
     'T_dwell1_0',
     'T_dwell2_0',
     'dwell2_pitch_set',
     'dwell2_times',
     'dwell2_cool_temp')



Now lets look at the ``coolstats`` array. This array is similar to
``hot_constraints['minus_yz'].dwell1_stats`` in that it includes
statistical data that is useful for plotting. The ``perc10``,
``perc50``, ``perc90`` fields contain the 10th, 50th, and 90th
percentile values for the cooldown times for each pitch. The cooldown
time denoted by the ``perc10`` column represents the time it will take
the constraint to reach the cooled state 10% of the time; this could be
considered an exceedingly optimistic (short) cooldown time. The
``perc90`` column denotes the time that it will take the constraint to
reach the cooled state 90% of the time; this would be considered the
conservative (longer) cooldown time.

In cases where factors other than attitude also impact temperature,
there may be a very large variation in cooldown time. In some of these
cases the 90th percentile cooldown time may reach the maximum dwell time
for the set of cooldown simulations, while the 10th percentile time will
be substantially lower. This is a real and expected behavior that stems
from the impacts of other nodes and power states (such as from the
number of ACIS CCDs used). In the case of models impacted by ACIS CCDs,
remember that if the MSID/constraint needs to be cooled following a hot
observation, the subsequent observation will use a reduced number of
CCDs, or an HRC observation will be scheduled.
In[116]:
.. code:: python

    coolstats.dtype.names





.. parsed-literal::
    ('pitch', 'perc10', 'perc50', 'perc90')



The ``hotstats`` array contains relevant statistical information
regarding the hot dwells used to seed each set of cooldown dwells. The
``dwell1_duration`` column contains the duration for each hot dwell. The
``dwell1_duration_delta`` contains the fraction of heatup time relevant
to the cooldown dwell times. The data in these two columns will be
explained in further detail below.

There is a very important difference between this information and the
data present in ``hot_constraints['minus_yz'].dwell1_stats``. The
``dwell1_stats`` array contains generalized statistical data and does
not contain data specific to each heatup dwell. This specific
information is necessary to set up the startup conditions for each set
of cooldown simulations. The ``hotstats`` array contains dwell duration
data for each individual hot dwell simulation picked to seed each set of
cooldown simulations, in the ``pitch`` and ``dwell1_duration`` columns.
The ``dwell1_duration_delta`` column requires more explanation on how
the set of hot dwells are selected to seed each set of cooldown
simulations.

Both the ``coolstats`` and ``hotstats`` arrays were originally generated
to make plotting relevant data simpler, compared to navigating the
primary dataset ``cool_constraints``. I am rethinking this strategy and
may eventually combine all three outputs of ``pyger.calc_constraints2``.

Much of the code used for the following discussion will go beyond what
is needed for normal analyses.
In[117]:
.. code:: python

    hotstats.dtype.names





.. parsed-literal::
    ('pitch', 'dwell1_duration', 'dwell1_duration_delta')



Before explaining what ``dwell1_duration_delta`` represents lets look at
how the cooldown dwell data is calculated. First we should look at all
the heatup dwells included in the original
``hot_constraints['minus_yz'].dwells1`` datastructure. This contains the
results of each node in the ``minus_yz`` model for each hot dwell
simulation, regardless of whether or not the temperature reached the
specified limit(s). The plot below shows the temperature for
``tcylaft6`` for each simulation.
In[144]:
.. code:: python

    # Create figure
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    # Plot each hot dwell simulation
    for d in hot_constraints['minus_yz'].dwells1:
        minus_yz_times = d.times
        tcylaft6_temps_deg_F = d.Ts[3] * 9 / 5 + 32
        ax.plot(minus_yz_times, tcylaft6_temps_deg_F, color=[0.8,0.8,0.8])
        
    # Format the plot
    ax.grid()
    ax.plot(ax.get_xlim(), [105, 105], 'b', linewidth=2)
    ax.set_title('TCYLAFT6 Temperatures for All Hot Simulations')
    ax.set_ylabel('Temperature in Deg F')
    ax.set_xlabel('Dwell Time in Seconds')





.. parsed-literal::
    <matplotlib.text.Text at 0x3281c550>




.. image:: index_files/index_52_1.png


We are not interested in simulating the cooldown time from dwells that
don't reach the specified limit, so the dwells that do not reach the
limit are filtered out as shown below.
In[145]:
.. code:: python

    # Create figure.
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    # Plot each hot dwell simulation if the max hot dwell duration is less than 90% of 
    # the max dwell time.
    # Filter out all dwell durations that show zero time.
    for d in hot_constraints['minus_yz'].dwells1:
        minus_yz_times = d.times
        tcylaft6_temps_deg_F = d.Ts[3] * 9 / 5 + 32
        dwell_duration = d.duration
        max_dwell_sec = hot_constraints['minus_yz'].max_dwell_ksec * 1000
        if (dwell_duration < max_dwell_sec) & (dwell_duration > 0.0):
            ax.plot(minus_yz_times, tcylaft6_temps_deg_F, color=[0.8,0.8,0.8])
            
    # Format the plot.
    ax.grid()
    ax.plot(ax.get_xlim(), [105, 105], 'b', linewidth=2)
    ax.set_title('TCYLAFT6 Temperatures for All Hot Simulations Reaching Limit')
    ax.set_ylabel('Temperature in Deg F')
    ax.set_xlabel('Dwell Time in Seconds')





.. parsed-literal::
    <matplotlib.text.Text at 0x42a07510>




.. image:: index_files/index_54_1.png


Since we want each cooldown dwell to be paired with a best effort hot
dwell, so that maximum likely cooldown time is comparable to maximum
likely heatup time, the hot dwell simulations that reach the specified
limit with the lowest starting temperatures are used. Specifically, of
the dwells that reach the specified limit, those with starting
temperatures at or under the 20th percentile are selected to seed the
cooldown dwells.
In[146]:
.. code:: python

    def find_lowest_starting_temp_indices(hot_constraints, model, msid, perc=0.2):
        """ Determine indices into dwells1 for hot dwells with lowest starting temps.
    
        :param hot_constraints: Output from pyger.calc_constraints() (hot dwells)
        :param model: Model of interest
        :param msid: MSID of interest
        :param perc: Bottom ratio of starting temperatures
        """
        
        msid_ind = hot_constraints[model].msids.index(msid)
        good_sim = hot_constraints[model].dwells1['duration'] > 0
        hot_pitch_ind_all = np.flatnonzero(
                                (hot_constraints[model].dwells1['constraint_name']
                                 != 'none') & good_sim)
        start_temps = hot_constraints[model].dwells1['T0s'][hot_pitch_ind_all]
        start_temps_transpose = start_temps.transpose()
        sort_ind = start_temps_transpose[msid_ind].argsort() 
        ind = range(int(len(start_temps_transpose[0]) * perc))
        lowest_temps_ind = sort_ind[ind]
        
        return hot_pitch_ind_all[lowest_temps_ind]
    
    # Call function to find indices of hot dwells of interest.
    msid = 'tcylaft6'
    model = 'minus_yz'
    hot_pitch_ind = find_lowest_starting_temp_indices(hot_constraints, model, 
                                                      msid, perc=0.2)
    
    # Create figure.
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    # Plot each hot dwell simulation if the max hot dwell duration is less than 90% of
    # the max dwell time.
    # Filter out all dwell durations that show zero time.
    for d in hot_constraints['minus_yz'].dwells1:
        minus_yz_times = d.times
        tcylaft6_temps_deg_F = d.Ts[3] * 9 / 5 + 32
        dwell_duration = d.duration
        max_dwell_sec = hot_constraints['minus_yz'].max_dwell_ksec * 1000 * 0.9
        if (dwell_duration < max_dwell_sec) & (dwell_duration > 0.0):
            ax.plot(minus_yz_times, tcylaft6_temps_deg_F, color=[0.8,0.8,0.8])
    
    # Plot hot dwells of interest in black.
    for d in hot_constraints['minus_yz'].dwells1[hot_pitch_ind]:
        minus_yz_times = d.times
        tcylaft6_temps_deg_F = d.Ts[3] * 9 / 5 + 32
        ax.plot(minus_yz_times, tcylaft6_temps_deg_F, color=[0,0,0]) 
    
    # Format plot.
    ax.grid()
    ax.plot(ax.get_xlim(), [105, 105], 'b', linewidth=2)
    ax.set_title('TCYLAFT6 Temperatures Used to Seed Cooldown Simulations (Black)')
    ax.set_ylabel('Temperature in Deg F')
    ax.set_xlabel('Dwell Time in Seconds')





.. parsed-literal::
    <matplotlib.text.Text at 0x42c3f3d0>




.. image:: index_files/index_56_1.png


Next look at what an example heatup/hot dwell and a subsequent cooldown
dwell. You'll see that in this case the cooldown dwell does not reach
the starting temperature observed at the begining of the heatup dwell.
There are some cases where the heatup dwell starting temperature can
easily be achieved, but this will not always be the case.

First find a good example cooldown dwell to use for this discussion
In[134]:
.. code:: python

    cool_pitch_vals = cool_constraints.dwell2_pitch_set[0][:]
    inds = np.arange(len(cool_pitch_vals))[(cool_pitch_vals > 130 & 
                                            (cool_pitch_vals < 135))]
    example_ind = inds[1]

The progression of temperatures during each cooldown dwell is not saved,
so to plot these tempemperatures, they need to be recalculated. This is
not normally something you would need to do, however it is useful for
this discussion.The following code recalculates these temperatures.
In[147]:
.. code:: python

    example_ind = inds[7]
    
    # Time values for the cooldown dwell (uses initial hot dwell time so date is approx. 
    # correct)
    start = hot_constraints['minus_yz'].start
    max_dwell_ksec = 600 # value used for pyger.constraints2() input
    stop = DateTime(start.secs + max_dwell_ksec * 1000)
    norm_profile = np.linspace(0, 1, 100)**10 # nonlinear time resolution
    times = start.secs + norm_profile * (stop.secs - start.secs)
    cool_times = times - times[0] # relative times useful for plotting
    
    # Obtain the temperature at which temperatures in the cooldown dwell are considered
    # 'cooled'
    T_cooled = cool_constraints.dwell2_cool_temp[example_ind] * 9 / 5 + 32
    
    # Obtain one of the target pitch values in the set of cooldown dwells
    cool_pitch = cool_constraints.dwell2_pitch_set[0][example_ind]
    
    # Re-use the original _get_states1 function to assemble states to calculate the 
    # cooldown temperatures
    states2 = hot_constraints['minus_yz']._get_states1(start, stop, cool_pitch)
    
    # Obtain the initial model temperatures for the example cooldown dwell
    T_dwell2_0 = cool_constraints.T_dwell2_0[example_ind]
    
    # Re-calculate the temperatures during the cooldown dwell example
    Ts_cool = hot_constraints['minus_yz']._calc_model(states2, times, T_dwell2_0)
    
    # Save information for the example hot seed dwell
    hot_ind = cool_constraints.dwell1_ind[example_ind]
    
    # Obtain temperatures during hot seed dwell
    hot_temps = hot_constraints['minus_yz'].dwells1[hot_ind].Ts[3]
    
    # Find indexes in hot_temps where temperatures are <= limit
    ind_under_limit = hot_temps <= (105 - 32) * 5 / 9
    
    # Add limit temperature to hot_temps to make plot flow better
    hot_temps = list(hot_temps[ind_under_limit] * 9 / 5 + 32)
    hot_temps.append(105)
    hot_temps = np.array(hot_temps)
    
    # Add duration to match temp. limit added above
    duration = hot_constraints['minus_yz'].dwells1[hot_ind].duration
    hot_times = hot_constraints['minus_yz'].dwells1[hot_ind].times - 
                hot_constraints['minus_yz'].dwells1[hot_ind].times[0]
    hot_times = list(hot_times[ind_under_limit])
    hot_times.append(duration)
    
    # Make hot_times start negative and move towards zero for plotting purposes
    hot_times = np.array(hot_times) - hot_times[-1]
    
    # Define x tick mark locations
    xtick = np.concatenate((np.linspace(hot_times[0], hot_times[-1], 6), 
                            np.linspace(cool_times[0], cool_times[-1], 8)), axis=1)
    
    # Define figure and plot relevant information
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    ax.plot(cool_times, Ts_cool[3]*9/5+32)
    ax.plot(hot_times, hot_temps,'r')
    ax.plot(ax.get_xlim(), [105, 105], 'r--', linewidth=2)
    ax.plot(ax.get_xlim(), [T_cooled, T_cooled], 'b--', linewidth=2)
    ax.plot(ax.get_xlim(), [hot_temps[0], hot_temps[0]], 'k--', linewidth=2)
    
    ax.set_xticks(xtick)
    ax.set_xlim((xtick[0], xtick[-1]))
    ax.grid()
    
    ax.text(40000, 105+.4, 'Constraint Limit', fontsize=12)
    ax.text(40000, T_cooled+0.4, 'Cooled Temperature Threshold', fontsize=12)
    ax.text(-225000, hot_temps[0]+0.4, 'Heatup Dwell Starting Temperature', fontsize=12)
    
    ax.set_title('Example Heatup and Cooldown Dwell')
    ax.set_ylabel('Temperature in Deg F')
    ax.set_xlabel('Time in Seconds')





.. parsed-literal::
    <matplotlib.text.Text at 0x442a1550>




.. image:: index_files/index_60_1.png


If you are running this notebook independantly you will have a different
random set of hot dwells and you may notice that the final temperature
is above the 'Cooled Temperature Threshold' or below the 'Heatup Dwell
Starting Temperature', however it is important to note that the latter
case will not always be common. For some constraints, the minimum
temperature will only be achievable at a cool attitude within a narrow
pitch range, or in the case of the ACIS related constraints, the
starting hot dwell starting temperatures as chosen above may not be
achievable. This latter case is possible when the hot dwell started
following an HRC observation, then 5 or 6 chips are used for both the
heatup and cooldown dwell.

Adding a pad in the form of the 'Cooled Temperature Threshold' shown
above expands the range that will reveal pitch ranges that will cool the
constraint. The downside of using such a pad is that it changes the
temperature range the reported cooldown time is relative to, so the
cooldown time is relative to a different temperature change compared to
the heatup dwell. Since we want the heatup time and the cooldown time to
be relative to the same change in temperature, the time between the
start of the heatup dwell and the point where the heatup dwell reaches
the 'Cooled Temperature Threshold' is removed from the heatup dwell time
and saved as the ``dwell1_duration_delta`` column in both the
``hotstats`` and ``cool_constraints`` arrays. The ratio of heatup time
between the values reported by ``dwell1_duration_delta`` and the
statistical values reported in the ``coolstats`` arrays are the values
that should be used to determine the ratio of required heating to
cooling time, such as when assembling the LTS.

Next a standard cooldown plot is generated using the built in method
``pyger.plot_cooldown``. The inputs to this function should be self
explanitory. For further details of this plot, please read the function
documentation; you can find the function in ``pyger.py``.

-  The black line represents the traditional best-case pyger maximum
   dwell capabiltiy in two degree pitch bins.
-  The orange datapoints represent the heatup time from the ``T_cooled``
   threshold to the limit (not from the starting temperature to the
   limit) for each hot dwell used to seed a set of cooldown simulations.
-  The orange line is the mean through these data.
-  The blue datapoints represent the amount of time required to reach
   the ``T_cooled`` threshold for each cooldown simulation.
-  The blue lines represent the ``perc10``, ``perc50``, and ``perc90``
   cooldown times that are comparable to the hot dwell data depicted in
   orange. The ``perc10`` line shows the amount of time required to cool
   this constraint 10% of the time (reports lowest time) where the
   \`perc90' shows the amount of time required to cool this constraint
   90% of the time (reports greatest time values)

You will see that the black line reaches a maximum value near 400 Ks,
where the cooldown data reaches 600 Ks, this is a result of using
different maximum dwell times in the heatup and cooldown simulation
runs. When creating your own plots you will probably want to customize
this plot to remove this discrepancy, however it is left here for
illustrative purposes.

Finally, to provide an example of how this plot should be used, consider
a hot dwell at a pitch of 70 degrees and a cool dwell at a pitch of 150
degrees. To understand the ratio of required heating to cooling time,
use the maximum dwell time at 70 degrees pitch as denoted by the orange
line and compare it to the cooling time at a pitch of 150 degrees using
the 90th percentile blue line (top line). Using the data available at
the time this example was generated, the max hot time of approximately
182 Ks would require a cooling time of 120 Ks to reach the original hot
dwell starting temperature with 90% confidence.
In[141]:
.. code:: python

    pyger.plot_cooldown(hot_constraints,
                        cool_constraints,
                        coolstats,
                        hotstats,
                        'minus_yz',
                        'tcylaft6',
                        '105.0',
                        'demo_tcylaft6.png')




.. image:: index_files/index_62_0.png

