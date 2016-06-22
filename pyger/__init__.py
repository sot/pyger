import warnings
warnings.simplefilter('ignore')

from .pyger import (calc_constraints, calc_constraints2, __version__, CtoF, ConstraintModel, 
                    ConstraintPSMC, ConstraintDPA, ConstraintDEA, ConstraintTank,
                    ConstraintAca, ConstraintFwdblkhd, ConstraintACISFP, save_pyger_pickle,
                    load_pyger_pickle, read_cases, run_pyger_cases, PostPyger)

from .make_sim_inputs import make_sim_inputs


def get_options(args=None):
    import argparse
    parser = argparse.ArgumentParser(prog='pyger')
    subparsers = parser.add_subparsers(title='Required sub-command (use exactly one)',
                                       help='Use "pyger <sub_command> -h" for help on sub-commands')

    # create the parser for the "make" command to make sim inputs
    parser_make = subparsers.add_parser('make', help='Make simulation inputs')
    parser_make.add_argument("--sim-file",
                             default='sim_inputs.pkl',
                             help="Output file (default = sim_inputs.pkl)")
    parser_make.add_argument("--start",
                             help="Sim input start time (default = Now - 1 year)")
    parser_make.add_argument("--stop",
                             help="Sim input stop time (default = Now - 10 days)")
    parser_make.add_argument("--n-days",
                             type=float,
                             help="Number of days to propagate prior to perigee exit (default = 3)")
    parser_make.add_argument("--min-dwell-sec",
                             type=float,
                             help="Minimum number of seconds to use as a propagation ending dwell (default = 1000)")
    parser_make.add_argument("--max-dwell-num",
                             type=float,
                             help="Maximum number of sample dwells to save (default = 300)")
    parser_make.set_defaults(func=make_sim_inputs)

    # create the parser for the "sim" command to run a constraint simulation
    parser_sim = subparsers.add_parser('sim', help='Simulate constraints')
    parser_sim.add_argument("--start",
                            help="Simulation start time")
    parser_sim.add_argument("--n-sim",
                            type=int,
                            help="Number of simulation points")
    parser_sim.add_argument("--dt",
                            type=float,
                            help="Delta time for sims (sec)")
    parser_sim.add_argument("--max-tcylaft6",
                            type=float,
                            help="TCYLAFT6 planning limit (degF)")
    parser_sim.add_argument("--max-1pdeaat",
                            type=float,
                            help="1PDEAAT planning limit (degC)")
    parser_sim.add_argument("--max-1dpamzt",
                            type=float,
                            help="1DPAMZT planning limit (degC)")
    parser_sim.add_argument("--max-pftank2t",
                            type=float,
                            help="PFTANK2T planning limit (degF)")
    parser_sim.add_argument("--max-4rt700t",
                            type=float,
                            help="4RT700T planning limit (degF)")
    parser_sim.add_argument("--max-fptemp-11",
                            type=float,
                            help="FPTEMP_11 limit (degC)")
    parser_sim.add_argument("--n-ccd",
                            type=int,
                            help="Number of ACIS CCDs")
    parser_sim.add_argument("--max-aacccdpt",
                            type=float,
                            help="ACA CCD planning limit (degC)")
    parser_sim.add_argument("--max-dwell-ksec",
                            type=float,
                            help="Max allowed dwell time (ksec)")
    parser_sim.add_argument("--sim-file",
                            default="sim_inputs.pkl",
                            help="Simulation inputs pickle file")
    parser_sim.set_defaults(func=calc_constraints)

    return parser.parse_args(args)


def main(args=None):
    import inspect
    opt = get_options(args)
    func_args = inspect.getargspec(opt.func)[0]
    opt_args = dict((k, v) for k, v in opt.__dict__.items()
                    if (k in func_args and v is not None))
    if opt.func == calc_constraints:
        import matplotlib
        matplotlib.use('Agg')
        constraints = calc_constraints(**opt_args)
    elif opt.func == make_sim_inputs:
        make_sim_inputs(**opt_args)

if __name__ == '__main__':
    main()
