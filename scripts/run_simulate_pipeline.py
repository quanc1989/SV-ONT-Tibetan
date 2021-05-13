# %%
import math
import os
import random
import sys

import msprime
import numpy as np
import tskit.vcf as vcf

import dadi
from Three_Population_Pipeline.Models_3D import split_sym_mig_size, split_symmig_all, split_symmig_all_grow, \
    sim_split_sym_mig_all
from my_models import mine_split_symmig_all, mine_split_symmig_all_anc

chrom_simulate = sys.argv[1]
label_simulate = sys.argv[2]
path_fs = sys.argv[3]
path_simulate = sys.argv[4]
model = sys.argv[5]
prefix = sys.argv[6]
method_random = sys.argv[7]
path_genetic_map = sys.argv[8]

assert chrom_simulate
assert method_random == 'uniform' or method_random == 'normal'

# path for genetic maps
length_window = 50000
# human mutation rate per site per generation
mu_custom = 1.5e-8
# human mutation rate per site per year
# mu = 0.5e-9

pts = [80, 90, 100]

# create an extrapolating function
if model == 'split_sym_mig_size':

    params_opt = [4.0624,0.3086,1.3716,0.4366,
                  0.0115,0.2576,1.0511,0.0009,
                  0.0218,0.0046,0.0262,0.0187,
                  0.0175,0.0364,0.0118,0.0053,
                  0.0023]
    theta_opt = 0

    uncerts = [0,0,0,0,
               0,0,0,0,
               0,0,0,0,
               0,0]

    GIM = []
    func_exec = dadi.Numerics.make_extrap_func(split_sym_mig_size)

elif model == 'split_symmig_all':

    if path_fs.__contains__('TIB_HANN_HANS'):

        params_opt = [3.5945, 4.7645, 4.9922, 0.3422, 0.004,
                      0.0212, 0.0991, 0.001, 0.0294, 0.0029]

        theta_opt = 96568.48

        uncerts = [1.53726012e-03, 9.38184756e-04, 3.66833417e-03, 3.38595564e-03, 5.09161035e-04,
                   4.24952519e-03, 7.51681635e-04, 2.38675729e-03, 1.24791601e-04, 4.10297270e-04,
                   1.28190791e+03]

        GIM_array = np.array([[3.29735533e+15,-2.23272909e+15,-7.59969166e+14, 6.79741119e+13,
                               2.26172867e+15,-1.03689428e+15,-2.66533767e+15, 3.11470816e+15,
                              -2.14014303e+15, 2.74059264e+15, -5.85046984e+07],
                             [-2.23272909e+15, 3.04932722e+15,  3.51423857e+15,  2.27782358e+15,
                              -1.43465778e+15, 2.52541035e+15,  3.84374887e+15, -1.04245093e+15,
                               4.13052141e+15, -3.58976745e+14,  5.19869627e+07],
                             [-7.59969166e+14, 3.51423857e+15,  9.48739602e+15,  9.68267640e+15,
                               2.83260997e+14, 5.86626127e+15, 8.75969445e+15, 5.87797954e+15,
                               6.06841173e+15, 3.70292498e+15,  4.70721121e+07],
                             [ 6.79741119e+13, 2.27782358e+15,  9.68267640e+15,  1.18126731e+16,
                               1.22531643e+15, 5.28701232e+15,  9.56464900e+15,  8.32936705e+15,
                               4.31165425e+15, 5.41226812e+15,  3.12163604e+07],
                             [ 2.26172867e+15, -1.43465778e+15,  2.83260997e+14,  1.22531643e+15,
                               1.71932565e+15, -3.56924219e+14, -9.12641139e+14,  2.95300754e+15,
                              -1.29688226e+15, 2.47707579e+15, -3.76306817e+07],
                             [-1.03689428e+15, 2.52541035e+15,  5.86626127e+15,  5.28701232e+15,
                              -3.56924219e+14, 4.29562662e+15,  5.49041732e+15,  3.17733703e+15,
                               4.27019955e+15, 7.10646915e+14,  3.89260013e+07],
                             [-2.66533767e+15, 3.84374887e+15,  8.75969445e+15,  9.56464900e+15,
                              -9.12641139e+14, 5.49041732e+15,  1.00626070e+16,  4.33072871e+15,
                               5.59746803e+15, 1.94042753e+15,  7.48177000e+07],
                             [ 3.11470816e+15, -1.04245093e+15,  5.87797954e+15,  8.32936705e+15,
                               2.95300754e+15,  3.17733703e+15,  4.33072871e+15,  9.74481877e+15,
                               3.48081724e+14,  5.11123755e+15, -3.45122480e+07],
                             [-2.14014303e+15,  4.13052141e+15,  6.06841173e+15,  4.31165425e+15,
                              -1.29688226e+15,  4.27019955e+15,  5.59746803e+15,  3.48081724e+14,
                               6.19281593e+15,  5.78571753e+14,  6.05143889e+07],
                             [ 2.74059264e+15, -3.58976745e+14,  3.70292498e+15,  5.41226812e+15,
                               2.47707579e+15,  7.10646915e+14,  1.94042753e+15,  5.11123755e+15,
                               5.78571753e+14,  6.25613088e+15, -3.30463391e+07],
                             [-5.85046984e+07,  5.19869627e+07,  4.70721121e+07,  3.12163604e+07,
                              -3.76306817e+07,  3.89260013e+07,  7.48177000e+07, -3.45122480e+07,
                               6.05143889e+07, -3.30463391e+07,  1.16448235e+00]])

        num_snv_reserved = 876639
        length_noncoding = 470094265
        num_snv_noncoding = 2565119

        ns = [2*80, 2*38, 2*46]

    elif path_fs.__contains__('YRI_TIB_HAN'):
        params_opt = [2.3118, 0.7128, 2.7737, 3.4771, 0.0166,
                      0.0983, 0.0072, 0.001, 0.0906, 0.0317]

        theta_opt = 130826.28

        uncerts = [8.62333129e-04,2.06568715e-03,2.11284160e-03,4.17208382e-04,1.41859306e-03,
                   1.39151011e-03,8.95924340e-04,4.06843195e-04,1.79688241e-03,2.26039772e-03,
                   3.01458991e+03]

        GIM_array = np.array([[2.49439744e+14, -1.59804356e+14, -2.98760675e+14, -1.50603121e+15,
                              -4.96010572e+14,  1.20095239e+15, -1.18871335e+14, -6.77165955e+14,
                               4.92925820e+14, -1.78962146e+15, -2.70731167e+06],
                             [-1.59804356e+14,  4.31117991e+14,  8.81010798e+14,  2.76757415e+15,
                               8.16149929e+14, -2.24828486e+15,  5.42300642e+14,  1.54388787e+15,
                              -3.90301441e+14,  3.46990137e+15,  6.21487208e+06],
                             [-2.98760675e+14,  8.81010798e+14,  1.97754257e+15,  6.35268586e+15,
                               1.84531851e+15, -5.10999399e+15,  1.16708443e+15,  3.57294296e+15,
                              -1.05364705e+15,  7.88424589e+15,  1.44155661e+07],
                             [-1.50603121e+15,  2.76757415e+15,  6.35268586e+15,  2.32242745e+16,
                               6.82085020e+15, -1.87456525e+16,  3.65015964e+15,  1.23650467e+16,
                              -4.86865426e+15,  2.81768046e+16,  5.10373396e+07],
                             [-4.96010572e+14,  8.16149929e+14,  1.84531851e+15,  6.82085020e+15,
                               2.03621074e+15, -5.47104750e+15,  1.02104818e+15,  3.60054573e+15,
                              -1.51718063e+15,  8.28114226e+15,  1.46775631e+07],
                             [ 1.20095239e+15, -2.24828486e+15, -5.10999399e+15, -1.87456525e+16,
                              -5.47104750e+15,  1.52362054e+16, -3.02807545e+15, -9.95249572e+15,
                               3.76890256e+15, -2.27206938e+16, -4.15257641e+07],
                             [-1.18871335e+14,  5.42300642e+14,  1.16708443e+15,  3.65015964e+15,
                               1.02104818e+15, -3.02807545e+15,  7.85410318e+14,  2.06769743e+15,
                              -3.98955619e+14,  4.53757766e+15,  8.70555197e+06],
                             [-6.77165955e+14,  1.54388787e+15,  3.57294296e+15,  1.23650467e+16,
                               3.60054573e+15, -9.95249572e+15,  2.06769743e+15,  6.75823512e+15,
                              -2.39641532e+15,  1.51246786e+16,  2.76447268e+07],
                             [ 4.92925820e+14, -3.90301441e+14, -1.05364705e+15, -4.86865426e+15,
                              -1.51718063e+15,  3.76890256e+15, -3.98955619e+14, -2.39641532e+15,
                               1.66701455e+15, -5.71930530e+15, -9.53586894e+06],
                             [-1.78962146e+15,  3.46990137e+15, 7.88424589e+15,  2.81768046e+16,
                               8.28114226e+15, -2.27206938e+16,  4.53757766e+15,  1.51246786e+16,
                              -5.71930530e+15,  3.43378561e+16,  6.20721267e+07],
                             [-2.70731167e+06,  6.21487208e+06,  1.44155661e+07,  5.10373396e+07,
                               1.46775631e+07, -4.15257641e+07,  8.70555197e+06,  2.76447268e+07,
                              -9.53586894e+06,  6.20721267e+07,  1.15501413e-01]])

        num_snv_reserved = 1702937
        length_noncoding = 470094265
        num_snv_noncoding = 4587735

        ns = [2*108, 2*80, 2*84]

    L = num_snv_reserved * (length_noncoding / num_snv_noncoding)
    func_exec = dadi.Numerics.make_extrap_func(split_symmig_all)

elif model == 'split_symmig_all_grow':

    params_opt = [4.0624,0.3086,1.3716,0.4366,
                  0.0115,0.2576,1.0511,0.0009,
                  0.0218,0.0046,0.0262,0.0187,
                  0.0175,0.0364,0.0118,0.0053,
                  0.0023]
    theta_opt = 0

    uncerts = [0,0,0,0,
               0,0,0,0,
               0,0,0,0,
               0,0]
    GIM = []
    func_exec = dadi.Numerics.make_extrap_func(split_symmig_all_grow)

elif model == 'sim_split_sym_mig_all':

    params_opt = [4.0624,0.3086,1.3716,0.4366,
                  0.0115,0.2576,1.0511,0.0009,
                  0.0218,0.0046,0.0262,0.0187,
                  0.0175,0.0364,0.0118,0.0053,
                  0.0023]
    theta_opt = 0

    uncerts = [0,0,0,0,
               0,0,0,0,
               0,0,0,0,
               0,0]
    GIM = []
    func_exec = dadi.Numerics.make_extrap_func(sim_split_sym_mig_all)

prefix = os.path.join(prefix, model+'.simulation', chrom_simulate)

if not os.path.exists(prefix):
    os.makedirs(prefix)

fs = dadi.Spectrum.from_file(path_fs)

# population size and grid points
popid = fs.pop_ids

# input params
print('optimized model: ', model)
print('populatin clusters: ', popid)
print('populatin projections: ', fs.sample_sizes)
print('populatin size to simulate: ', ns)
print('optimized theta: ', theta_opt)
print('optimized parameters : ', params_opt)
print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

#read fs for each chromosome and calculate theta for each windows in order
# example:
# theta_all_windows = {
#     'chr22': [1.332, 2.664]
# }
# theta_max = {
#     'chr22': 2.664
# }


dict_fs = {}

# store dict fs per window
for chrom in os.listdir(path_simulate):

    if chrom == chrom_simulate:

        dict_fs[chrom] = {}

        # read fs for each window
        # index_fs = 0
        for window in os.listdir(os.path.join(path_simulate, chrom)):

            # if index_fs % 1000 == 0:
            #     print(chrom, index_fs)
            # index_fs += 1
            # print(path_simulate, chrom, window)
            fs_window = dadi.Spectrum.from_file(
                os.path.join(path_simulate, chrom, window)
            )

            pos = int(window.split('.')[0].split('-')[1])
            dict_fs[chrom][pos] = fs_window


def generate_random_parameters(list_params, init_theta, list_uncerts, raw_sfs, raw_pts, dict_sfs, GIM=None):

    if GIM is None:
        # method1: uniform draw from confidence interval
        index_param = 0
        params_opt_random = [x for x in list_params]
        while index_param < len(params_opt_random):
            # params_opt_random[index_param] = np.random.normal(params_opt[index_param], uncerts[index_param])
            params_opt_random[index_param] = np.random.uniform(list_params[index_param]-1.96*list_uncerts[index_param],
                                                               list_params[index_param]+1.96*list_uncerts[index_param])
            index_param += 1
        theta_random = np.random.uniform(init_theta-1.96*list_uncerts[-1],
                                         init_theta+1.96*list_uncerts[-1])
    else:
        # method2: multivariate-normal distribution, need covariance matrix
        tmplist = [x for x in list_params]
        tmplist.append(init_theta)

        cov_matrix = np.linalg.inv(GIM)
        params_opt_random = np.random.multivariate_normal(tmplist, cov_matrix)
        theta_random = params_opt_random[-1]
        params_opt_random = params_opt_random[:-1]

    Ne_random = theta_random / (4 * mu_custom * L)
    sim_model = func_exec(params_opt_random, raw_sfs.sample_sizes, raw_pts)

    theta_all_windows_segments = {}
    for chrom in dict_sfs:
        theta_all_windows_segments[chrom] = {}
        for pos in dict_sfs[chrom]:
            fs_window = dict_sfs[chrom][pos]
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs_window)
            theta_all_windows_segments[chrom][pos] = theta

    return params_opt_random, Ne_random, theta_all_windows_segments

# simulate the model with the optimized parameters
# for index_simulate in range(num_simulations):


if os.path.exists(os.path.join(prefix, "variants.%s.%s.filter.vcf"%(chrom, label_simulate))):
    print('simulation already exist')
else:
    # use GIM uncerts to random generate parameters in normal distribution for each simulation
    if method_random == 'uniform':
        params_opt_random, Ne_random, theta_all_windows_segments = generate_random_parameters(params_opt, theta_opt, uncerts, fs, pts, dict_fs)
    else:
        params_opt_random, Ne_random, theta_all_windows_segments = generate_random_parameters(params_opt, theta_opt, uncerts, fs, pts, dict_fs, GIM_array)

    print('simulation label ' + label_simulate)
    print('random parameters:' + ','.join([str(x) for x in params_opt_random]))
    print('original parameters:' + ','.join([str(x) for x in params_opt]))
    print('random population size: %f' % Ne_random)

    # do simulation for each chromosome using its maximum theta
    for chrom in theta_all_windows_segments:
        # calculate max theta
        theta_all_windows = [theta_all_windows_segments[chrom][x] for x in theta_all_windows_segments[chrom]]
        theta = max(theta_all_windows)
        # diploid effective populations calculated by theta/ (4*mu*L)
        # Ne = theta / (4 * mu * length_window)
        mu = theta / (4 * Ne_random * length_window)

        print('max segment theta: %f' % theta)
        print('estimated segment mu: %e' % mu)

        # Read in the recombination map using the read_hapmap method in msprime
        infile_genetic_map = os.path.join(path_genetic_map, 'genetic_map_GRCh37_%s.txt' % chrom)
        recomb_map = msprime.RecombinationMap.read_hapmap(infile_genetic_map)

        # define demographic model for msprime
        if model == 'split_symmig_all':
            if path_fs.__contains__('TIB_HANN_HANS'):
                population_configurations, migration_matrix, demographic_events = mine_split_symmig_all(Ne_random,
                                                                                                     params_opt_random,
                                                                                                     ns)
            elif path_fs.__contains__('YRI_TIB_HAN'):
                population_configurations, migration_matrix, demographic_events = mine_split_symmig_all_anc(Ne_random,
                                                                                                        params_opt_random,
                                                                                                        ns)
        else:
            print('error -> not split_symmig_all')

        # do simulate in a whole chromosome
        tree_sequence = msprime.simulate(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events,
            Ne=Ne_random,
            mutation_rate=mu,
            recombination_map=recomb_map,
        )

        # save simulated vairants as VCF
        with open(os.path.join(prefix, "variants.%s.%s.vcf" % (chrom, label_simulate)), "w") as vcf_file:
            writer = vcf.VcfWriter(tree_sequence, ploidy=2, contig_id=chrom)
            writer.write(vcf_file)

        cache_label = 0
        cache_label_pos = []
        res_pos = []
        # drop a proportion of (1- theta/theta_max) of variants for each window
        with open(os.path.join(prefix, "variants.%s.%s.theta" % (chrom, label_simulate)), "w") as theta_file:

            theta_file.write('%s\tMax\t%f\n' % (chrom, theta))

            for variant in tree_sequence.variants():
                label = int(variant.site.position/length_window)
                if label != cache_label:
                    if cache_label in theta_all_windows_segments[chrom]:
                        if len(cache_label_pos) > 0:
                            theta_segment = theta_all_windows_segments[chrom][cache_label]
                            maxnum_variants = math.ceil(len(cache_label_pos) * (1.0 * theta_segment) / theta)

                            # print('segment theta %f' %theta_segment)
                            # print('maxnum variants %f' % maxnum_variants)
                            # print('cache_label_pos %d' % len(cache_label_pos))

                            if maxnum_variants > len(cache_label_pos):
                                maxnum_variants = len(cache_label_pos)
                            res_pos.extend(random.sample(cache_label_pos, k=maxnum_variants))
                            theta_file.write('%s\t%d\t%d\t%f\t%d\n' % (chrom, cache_label, len(cache_label_pos), theta_segment, maxnum_variants))
                        else:
                            theta_file.write('%s\t%d\t%d\t%f\t%d\n' % (chrom, cache_label, len(cache_label_pos), theta_segment, maxnum_variants))
                    else:
                        theta_file.write('%s\t%d\t%d\t0\t0\n' % (chrom, cache_label, len(cache_label_pos)))
                    cache_label = label
                    cache_label_pos = []

                else:
                    cache_label_pos.append(variant.site.position)

            if len(cache_label_pos) > 0:
                if cache_label in theta_all_windows_segments[chrom]:
                    theta_segment = theta_all_windows_segments[chrom][cache_label]
                    maxnum_variants = math.ceil(len(cache_label_pos) * (1.0 * theta_segment) / theta)
                    res_pos.extend(random.sample(cache_label_pos, k=maxnum_variants))
                    theta_file.write(
                        '%s\t%d\t%d\t%f\t%d\n' % (chrom, cache_label, len(cache_label_pos), theta_segment, maxnum_variants))
                else:
                    theta_file.write('%s\t%d\t%d\t0\t0\n' % (chrom, cache_label, len(cache_label_pos)))

        with open(os.path.join(prefix, "variants.%s.%s.pos" % (chrom, label_simulate)), "w") as pos_file:
            for pos in res_pos:
                pos_file.write('%s\t%d\n' % (chrom, pos))

        command_vcftools = 'vcftools --vcf ' + os.path.join(prefix, "variants.%s.%s.vcf"%(chrom, label_simulate)) + \
                           ' --recode --positions ' + os.path.join(prefix, "variants.%s.%s.pos"%(chrom, label_simulate)) + \
                           ' --stdout > ' + os.path.join(prefix, "variants.%s.%s.filter.vcf"%(chrom, label_simulate))
        print(command_vcftools)
        os.system(command_vcftools)
        command_rm = 'rm ' + os.path.join(prefix, "variants.%s.%s.pos"%(chrom, label_simulate))
        print(command_rm)
        os.system(command_rm)
        os.system('bgzip ' + os.path.join(prefix, "variants.%s.%s.vcf"%(chrom, label_simulate)))
        os.system('bgzip ' + os.path.join(prefix, "variants.%s.%s.filter.vcf"%(chrom, label_simulate)))
