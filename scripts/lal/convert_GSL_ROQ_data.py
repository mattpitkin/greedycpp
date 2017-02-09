#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""Script to convert GSL binary data output by greedycpp and the eim to numpy or hdf5.

This script reads data for the reduced basis Basis.bin', 'ApproxErrors.txt',
'quad_rule.txt', the empirical interpolant data stored either in
'invV.bin' and 'EIM_indices.txt' or 'Aminv.bin' and 'EIM_X.txt'.

It sorts the EIM frequencies and B matrix so it plays nice with assumptions 
made in the LAL FrequencySequence interface to the Fourier domain waveforms.

It saves data in numpy or hdf5 (preferred) format. The advantage of hdf5
is that you get everything in one file including metadata.

The training set file is read and bounds for physical parameters are saved as metadata.
I assume that the basis was created by sampling in chirp mass and symmetric mass-ratio,
among other parameters, and save these bounds.

Numpy style output creates these files::
    'EIM_F_sorted.npy'
    'EIM_dfs_sorted.npy'
    'B_sorted.npy'

Example usage. We have a bunch of ROQ rules for the linear/quadratic parts
of a likelihood expression, and in different mass/time bands.

Write a script called GatherROQ.sh to collect the data (written in gsl binary)
and place it into a common folder called PhenomPv2ROQ

>>> cat GatherROQ.sh

#!/bin/bash

# 4s #
DIR='/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_PlusCross_parts_4s/it1/'
./convert_GSL_ROQ_data.py -d ${DIR} -s numpy -c ${DIR}run_settings.cfg -p ./PhenomPv2ROQ/4s/linear
DIR='/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_Quad_parts_4s/'
./convert_GSL_ROQ_data.py -d ${DIR} -s numpy -c ${DIR}run_settings.cfg -p ./PhenomPv2ROQ/4s/quadratic

# 8s #
DIR='/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_PlusCross_parts_8s/it1/'
./convert_GSL_ROQ_data.py -d ${DIR} -s numpy -c ${DIR}run_settings.cfg -p ./PhenomPv2ROQ/8s/linear
DIR='/panfs/ds08/sxs/sfield/greedy_phenomPv2/phenomp_Quad_parts_8s/'
./convert_GSL_ROQ_data.py -d ${DIR} -s numpy -c ${DIR}run_settings.cfg -p ./PhenomPv2ROQ/8s/quadratic

"""

__author__ = "Michael Puerrer"
__copyright__ = "Copyright 2015"
__email__ = "Michael.Puerrer@ligo.org"

import os.path
import numpy as np
import re, sys, os
from optparse import OptionParser

def etafun(q):
    return q/(1.0 + q)**2
def Mchirpfun(M, eta):
    return M*eta**(3.0/5.0)

def ReadComplexGSLBinaryMatrix(directory, filename, m, N):
    data = np.fromfile(os.path.join(directory, filename))
    # 'data' is an array of size 2 * m * N saved by gsl_matrix_complex_fwrite()
    # Entries are [Re(0,0), Im(0,0), Re(0,1), Im(0,1), ..., Re(m-1,N-1), Im(m-1,N-1)]
    A_re = np.reshape(data[::2], (m,N))
    A_im = np.reshape(data[1::2], (m,N))
    return A_re + 1.j*A_im

def FindTrainingRegion(TS_file, approximant):
    """ Deduce basis training region so that this information can be exported."""

    print 'Loading training set from file', TS_file
    TS = np.loadtxt(TS_file)
    num_pars = np.shape(TS)[1]
    if 'SEOBNRv2_ROM_DoubleSpin_HI' in approximant: # catch all the different polarization tags
        # The ROM interface used parameters ['m1', 'm2', 'chi1', 'chi2']
        m1 = TS.T[0]
        m2 = TS.T[1]
        etas = etafun(m1 / m2)
        Mcs = Mchirpfun(m1 + m2, etas)
        TS.T[0] = Mcs
        TS.T[1] = etas
        par_names = ['Mc', 'eta', 'chi1', 'chi2']
    elif approximant=='LackeyTidal2013_SEOBNRv2_ROM_HI_all_parts':
        # The model interface used parameters ['mBH', 'mNS', 'chiBH', 'Lambda']
        # The TS also uses these parameters, so we will stick to them to preserve the boundary.
        par_names = ['mBH', 'mNS', 'chiBH', 'Lambda']
    elif 'PhenomP' in approximant: # catch all the different polarizations and h^2 terms
        # The PhenomP interface used ['Mtot', 'eta', 'chi_eff', 'chip', 'thetaJ', 'alpha0']
        Mtots = TS.T[0]
        etas = TS.T[1]
        TS.T[0] = Mchirpfun(Mtots, etas)
        par_names = ['Mc', 'eta', 'chi_eff', 'chip', 'thetaJ', 'alpha0']

    return par_names, TS


def ParseConfigFile(cfg_file):

    model_pattern = re.compile('model_name[=\s]+"(\w+)"')
    for line in open(cfg_file):
        for match in re.finditer(model_pattern, line):
            approximant = match.groups()[0]
            break
    ts_pattern = re.compile('ts_file[=\s]+"(.*)"')
    for line in open(cfg_file):
        for match in re.finditer(ts_pattern, line):
            TS_file = match.groups()[0]
            break

    roq_partition_pattern = re.compile('roq_[0-9\s]')
    roq_partition = None
    for line in open(cfg_file):
        x = re.search(roq_partition_pattern, line)
        if x is not None:
            roq_partition = int(x.group().split('_')[1])
            break

    print "roq partition is %s"%roq_partition

    return approximant, TS_file, roq_partition

def ParseROQInfo(roq_info,outdir,roq_partition):
    """ Use roq_info.json file to deduce the correct params.dat file"""

    import json
    with open(roq_info) as fp:
        opts = json.load(fp)

    roq_partitions = opts["partitions"]
    ts_boundary_file = roq_partitions[roq_partition]["sample intervals"]

    x = np.genfromtxt(ts_boundary_file, names=True)

    par_names = '#'
    par_value = ''

    for param_name in x.dtype.names:
      par_names += param_name+'-min '+ param_name+'-max '
      par_min = x[param_name][0]
      par_max = x[param_name][1]
      par_value += str(par_min) + ' '+ str(par_max) + ' '

    par_names +='\n' 

    fp = open(os.path.join(outdir, 'params.dat'),'w')
    fp.write(par_names)
    fp.write(par_value)
    fp.close()

def Convert_GSL_EIM_Data(directory, hdf5_filename='ROQ_SEOBNRv2_ROM_LM_40_4096Hz.hdf5', 
      output_style='hdf5', cfg_file='run_40_4096Hz_pbnd_phase_corr_fine_local.cfg',
      outdir='.',roq_info=None):

    print(outdir)

    if directory != outdir and output_style == 'numpy':
        os.makedirs(outdir)

    # Read RB and EIM from directory, write B matrix for ROQ
    print 'Reading data from directory', directory
    ApproxErrors = np.loadtxt(os.path.join(directory, 'ApproxErrors.txt'))
    [freqs, dfs] = np.loadtxt(os.path.join(directory, 'quad_rule.txt')).T
    m = len(ApproxErrors) # = dim(RB)
    N = len(freqs) # number of frequency points
    print 'dim(RB) = ', m
    print 'dim(f_i) = ', N

    e = ReadComplexGSLBinaryMatrix(directory, 'Basis.bin', m, N)

    if os.path.isfile(os.path.join(directory, 'invV.bin')):
        print 'Reading EIM output by new O(n^3) eim code'
        Vinv = ReadComplexGSLBinaryMatrix(directory, 'invV.bin', m, m)
        EIM_idx = np.loadtxt(os.path.join(directory, 'EIM_indices.txt'), dtype=int)
    else:
        print 'Reading EIM output by old O(n^4) eim code'
        Vinv = ReadComplexGSLBinaryMatrix(directory, 'Aminv.bin', m, m)
        EIM_idx = np.loadtxt(os.path.join(directory, 'EIM_X.txt'), dtype=int)

    print 'Sorting EIM frequencies and B matrix'
    B_unsorted = np.dot(e.T, Vinv)
    EIM_F_unsorted = freqs[EIM_idx]
    EIM_dfs_unsorted = dfs[EIM_idx]

    # Note: EIM_F are unsorted, but the phase correction in SEOBNR ROMs and IMRPhenomP
    # only really works if the frequency sequence is sorted. If we don't sort, we'll 
    # potentially get a wrong coalescence time and thus significant error in the overlap.
    idx = np.argsort(EIM_F_unsorted)   # get sort keys
    EIM_F = EIM_F_unsorted[idx]        # sort Fs
    B = B_unsorted[:,idx]              # sort the B matrix
    EIM_dfs = EIM_dfs_unsorted[idx]

    # Parse greedcpp config file
    approximant, TS_file, roq_partition = ParseConfigFile(cfg_file)

    # Write training set bound metadata for all physical parameters
    par_names, TS = FindTrainingRegion(TS_file, approximant)

    # Parse roq_info.json file
    if roq_info is not None:
        ParseROQInfo(roq_info, outdir, roq_partition)
    else:
        print "No roq_info.json file to parse"
        

    if output_style=='numpy':
        print 'Saving B matrix and EIM in numpy format'
        print "model = %s"%approximant
        print "TS_file = %s"%TS_file
        np.save(os.path.join(outdir, 'EIM_F_sorted.npy'), EIM_F)
        np.save(os.path.join(outdir, 'EIM_dfs_sorted.npy'), EIM_dfs)
        np.save(os.path.join(outdir, 'B_sorted.npy'), B)
    elif output_style=='hdf5':
        if not hdf5_filename:
            hdf5_filename = 'ROQ' + '_' + approximant + '_' + str(np.min(freqs)) + '_' + str(np.max(freqs)) + 'Hz.hdf5'
        print 'Saving ROQ data in hdf5 format to', hdf5_filename
        import h5py
        fp = h5py.File(hdf5_filename, 'w')

        dset_B = fp.create_dataset('B', data=B)
        dset_B.attrs['Comment'] = np.string_('The sorted B matrix: B = e^T . inv(V). shape(B) = C^{N x m}')
        # Don't save e as it can be reconstructed from B and Vinv and is usually very large.
        #dset_e = fp.create_dataset('Reduced basis', data=e)
        #dset_e.attrs['Comment'] = np.string_('The reduced basis e. shape(e) = C^{m x N}')
        dset_Vinv = fp.create_dataset('V inverse', data=Vinv)
        dset_Vinv.attrs['Comment'] = np.string_('The inverse Vandermonde matrix. shape(inv(V)) = C^{m x m}')

        dset_ApproxErrors = fp.create_dataset('ApproxErrors', data=ApproxErrors)
        dset_EIM_dfs = fp.create_dataset('EIM_deltaF_j', data=EIM_dfs)
        dset_EIM_F = fp.create_dataset('EIM_F_j', data=EIM_F)
        dset_EIM_idx = fp.create_dataset('EIM_indices', data=idx)
        dset_grid_f_i = fp.create_dataset('f_i', data=freqs)
        dset_grid_df_i = fp.create_dataset('deltaf_i', data=dfs)

        print 'Assuming approximant', approximant
        fp.attrs['Approximant'] = np.string_(approximant)
        fp.attrs['f_min'] = np.min(freqs)
        fp.attrs['f_max'] = np.max(freqs)
        fp.attrs['dim(e)'] = m
        fp.attrs['dim(f_i)'] = N
        fp.attrs['Comment'] = np.string_(
            '''
            The frequencies EIM_F_j and frequency spacing EIM_deltaF_j of the empirical 
            interpolant have been sorted to be monotonically increasing as 
            required by LAL codes.
            The columns of the B matrix (= e^T . inv(V)) are sorted accordingly.
            
            We use the following abbreviations:
            m = dim(e)
            N = dim(f_i)
            
            e is the reduced basis.
            f_i are the frequencies of the quadrature rule used to build the ROQ
            and deltaf_i are the corresponding frequency spacings.
            '''
            )

        for i in range(len(par_names)):
            fp.attrs[par_names[i]+'_min'] = np.min(TS.T[i])
            fp.attrs[par_names[i]+'_max'] = np.max(TS.T[i])
        print 'Saved training set bounds for parameters', par_names 

        fp.close()
    else:
      print 'Unknown output format', output_style


if __name__ == "__main__":
    parser = OptionParser(description=__doc__)
    parser.add_option("-d", "--path-to-ROQ-data", dest="data_directory", type="string",
                  help="Path to directory containing GSL binary data for ROQ output by greedycpp and eim code.", 
                  metavar="data_directory")
    parser.add_option("-o", "--hdf5-filename", dest="hdf5_filename", type="string",
                  help="Optional name of the hdf5 output file. The name will be aut generated if not specified.", 
                  metavar="hdf5_filename")
    parser.add_option("-s", "--output-style", dest="output_style", type="string",
                  help="Output style: numpy or hdf5 (default).", metavar="output_style")
    parser.add_option("-c", "--cfg-file", dest="cfg_file", type="string",
                  help="Path to the greedycpp configuration file used to build the ROQ.", metavar="cfg_file")
    parser.add_option("-p", "--path-of-output", dest="path_of_output", type="string",
                  help="Directory path in which to place the output data.", 
                  metavar="data_directory")
    parser.add_option("-j", "--roq-info", dest="roq_info", type="string",
                  help="Path to the roq_info.json file, if it exists (this file is used in the roq pipeline). If provided, a params.dat fill will be created.", metavar="roq_info")

    parser.set_defaults(data_directory='.', hdf5_filename=None,
                        output_style='hdf5', cfg_file=None,path_of_output=None,
                        roq_info=None)

    (options, args) = parser.parse_args()

    if not options.cfg_file:
        print 'Please specify a greedcpp configuration file.'
        sys.exit(-1)

    # Read and convert data
    if options.path_of_output is None:
        outdir = options.data_directory
    else:
        outdir = options.path_of_output
    Convert_GSL_EIM_Data(options.data_directory, hdf5_filename=options.hdf5_filename, 
                         output_style=options.output_style, cfg_file=options.cfg_file,
                         outdir=outdir, roq_info=options.roq_info)

    print 'All done!\n'

