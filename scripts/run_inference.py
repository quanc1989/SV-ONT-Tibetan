import configargparse
import numpy

import Optimize_Functions
import dadi


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-s', '--sfs', required=True,
               help='sfs')
    parser.add('-p', '--prefix', required=True,
               help='save prefix')
    parser.add('-m', '--model', required=True,
               help='model name')
    parser.add('--pts', required=False, default='80,90,100',
               help='init pts')
    parser.add('--project', required=False,
               help='project to a small size')
    parser.add('--param', required=False,
               help='param')
    parser.add('--rounds', required=False, default=9, type=int,
               help='round number')
    parser.add('--unfolded', default=False, action='store_true',
               help='if unfolded')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    path_sfs = opt.sfs
    fs = dadi.Spectrum.from_file(path_sfs)
    if opt.project:
        fs = fs.project([int(x.strip()) for x in opt.project.split(',')])

    opt.folded = not opt.unfolded
    ns = fs.sample_sizes
    popid = fs.pop_ids
    print(fs.pop_ids)

    print("\n\n============================================================================")
    print("\nData for site frequency spectrum:\n")
    print("Projection: {}".format(ns))
    print("Sample sizes: {}".format(fs.sample_sizes))
    print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
    print("\n============================================================================\n")

    from Three_Population_Pipeline.Models_3D import split_nomig
    from Three_Population_Pipeline.Models_3D import split_symmig_all
    from Three_Population_Pipeline.Models_3D import split_symmig_all_init
    from Three_Population_Pipeline.Models_3D import split_symmig_all_grow
    from Three_Population_Pipeline.Models_3D import split_symmig_all_grow_init
    from Three_Population_Pipeline.Models_3D import split_asymmig_all
    from Three_Population_Pipeline.Models_3D import split_nomig_size
    from Three_Population_Pipeline.Models_3D import split_sym_mig_size
    from Three_Population_Pipeline.Models_3D import split_asym_mig_size
    from Three_Population_Pipeline.Models_3D import sim_split_no_mig
    from Three_Population_Pipeline.Models_3D import sim_split_no_mig_size
    from Three_Population_Pipeline.Models_3D import sim_split_sym_mig_all
    from Three_Population_Pipeline.Models_3D import sim_split_asym_mig_all
    from Three_Population_Pipeline.Models_3D import sim_split_refugia_sym_mig_all

    # create a prefix to label the output files
    prefix = opt.prefix

    # Remember the order for mandatory arguments as below
    # Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded)
    pts_l = [int(x.strip()) for x in opt.pts.split(',')]

    # reps = [10, 10, 20]
    # maxiters = [5, 5, 5]
    # folds = [3, 2, 1]

    if opt.param:
        inparam = [float(x) for x in opt.param.split(',')]
    else:
        inparam=None

    if opt.model == 'sim_split_no_mig':

        p_labels = 'nu1, nu2, nu3, T1'
        in_upper = [5, 5, 5, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.001]

        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="sim_split_no_mig",
                                            func=sim_split_no_mig,
                                            rounds=opt.rounds,
                                            param_number=4,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_nomig':
        p_labels = 'nu1, nuA, nu2, nu3, T1, T2'
        in_upper = [5, 5, 5, 5, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_nomig",
                                            func=split_nomig,
                                            rounds=opt.rounds,
                                            param_number=6,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)
    elif opt.model == 'split_nomig_size':
        p_labels = 'nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2, T3'
        in_upper = [5, 5, 5, 5, 5, 5, 5, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_nomig_size",
                                            func=split_nomig_size,
                                            rounds=opt.rounds,
                                            param_number=10,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_sym_mig_size':
        p_labels = 'nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, mA12, mB12, mB23, mB13, mC12, mC23, mC13, T1, T2, T3'
        in_upper = [5, 5, 5, 5, 5, 5, 5, 0.1, 0.1, 0.1,0.1, 0.1, 0.1,0.1, 0.1, 0.1,0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001,0.001, 0.001, 0.001,0.001, 0.001, 0.001,0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_sym_mig_size",
                                            func=split_sym_mig_size,
                                            rounds=opt.rounds,
                                            param_number=17,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_asym_mig_size':
        p_labels = 'nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, mA12, mA21, mB12, mB21, mB23, mB32, mB13, mB31, mC12, mC21, mC23, mC32, mC13, mC31, T1, T2, T3'
        in_upper = [5, 5, 5, 5, 5, 5, 5, 0.1, 0.1, 0.1,0.1, 0.1, 0.1,0.1, 0.1, 0.1,0.1, 0.1, 0.1,0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                    0.001, 0.001, 0.001,0.001, 0.001, 0.001,0.001,
                    0.001, 0.001,0.001, 0.001, 0.001,0.001, 0.001, 0.001,0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_asym_mig_size",
                                            func=split_asym_mig_size,
                                            rounds=opt.rounds,
                                            param_number=24,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_symmig_all':
        p_labels = 'nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2'
        in_upper = [5, 5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_symmig_all",
                                            func=split_symmig_all,
                                            rounds=opt.rounds,
                                            param_number=10,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_symmig_all_init':
        p_labels = 'nu0, nu1, nuA, nu2, nu3, mA, m1, m2, m3, T0, T1, T2'
        in_upper = [5, 5, 5, 5, 5,
                    0.1, 0.1, 0.1, 0.1, 0.1,
                    0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01,
                    0.001, 0.001, 0.001, 0.001, 0.0001,
                    0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_symmig_all_init",
                                            func=split_symmig_all_init,
                                            rounds=opt.rounds,
                                            param_number=12,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_symmig_all_grow':
        p_labels = 'nu1, nuA, nu2, nu3, nu1B, nu2B, nu3B, mA, m1, m2, m3, T1, T2'
        in_upper = [5, 5, 5, 5, 5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_symmig_all_grow",
                                            func=split_symmig_all_grow,
                                            rounds=opt.rounds,
                                            param_number=13,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_symmig_all_grow_init':
        p_labels = 'nu0, nu1, nuA, nu2, nu3, nu2B, nu3B, mA, m1, m2, m3, T0, T1, T2'
        in_upper = [5, 5, 5, 5, 5,
                    5, 5, 0.1, 0.1, 0.1,
                    0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01,
                    0.01, 0.01, 0.001, 0.001, 0.001,
                    0.001, 0.0001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_symmig_all_grow_init",
                                            func=split_symmig_all_grow_init,
                                            rounds=opt.rounds,
                                            param_number=14,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'split_asymmig_all':
        p_labels = 'nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m23, m32, m13, m31, T1, T2'
        in_upper = [5, 5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="split_asymmig_all",
                                            func=split_asymmig_all,
                                            rounds=opt.rounds,
                                            param_number=14,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'sim_split_sym_mig_all':
        p_labels = 'nu1, nu2, nu3, m1, m2, m3, T1'
        in_upper = [5, 5, 5, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="sim_split_sym_mig_all",
                                            func=sim_split_sym_mig_all,
                                            rounds=opt.rounds,
                                            param_number=7,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)

    elif opt.model == 'sim_split_refugia_sym_mig_all':
        p_labels = 'nu1, nu2, nu3, m1, m2, m3, T1, T2'
        in_upper = [5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="sim_split_refugia_sym_mig_all",
                                            func=sim_split_refugia_sym_mig_all,
                                            rounds=opt.rounds,
                                            param_number=8,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper
                                            )

    elif opt.model == 'sim_split_no_mig_size':
        p_labels = 'nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2'
        in_upper = [5, 5, 5, 5, 5, 5, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="sim_split_no_mig_size",
                                            func=sim_split_no_mig_size,
                                            rounds=opt.rounds,
                                            param_number=8,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper
                                            )

    elif opt.model == 'sim_split_asym_mig_all':
        p_labels = 'nu1, nu2, nu3, m12, m23, m31, m21, m32, m13, T1'
        in_upper = [5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        in_lower = [0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        Optimize_Functions.Optimize_Routine(fs=fs,
                                            pts=pts_l,
                                            outfile=prefix,
                                            model_name="sim_split_asym_mig_all",
                                            func=sim_split_asym_mig_all,
                                            rounds=opt.rounds,
                                            param_number=10,
                                            fs_folded=opt.folded,
                                            param_labels=p_labels,
                                            optimizer="log",
                                            # reps=reps,
                                            # maxiters=maxiters,
                                            # folds=folds,
                                            in_params=inparam,
                                            in_lower=in_lower,
                                            in_upper=in_upper)
