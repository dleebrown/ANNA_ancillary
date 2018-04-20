"""
v1.0 This implements the postprocessing routines - smoothing to a lower resolution, rotation, SN, radial velocity
(interpolation)
More or less built from specgen, so everything is basically the same. Reads in high resolution spectra output
by specgen, and then resamples from postprocessing parameter ranges and outputs specified number of files,
repeating input files if necessary.
Has the ability to output binary files with the format
[#px, #params, wavelengths, star1 params, star 1 flux, star2 params...] and remove *.vsn.spec files to save space

Inputs: *.spec files output by specgen, 'postproc.param' or parameter file
Outputs: *.vsn.spec files, binary output ('allspec')
"""

import subprocess
import random
import multiprocessing
from functools import partial
import time
import pandas as pd
import numpy as np
from array import array
import glob
import math
from scipy import interpolate as ip

# OKAY 05/25/16
# pass this function an open csv file, and it will write the starname and atmospheric parameters used
# just doing this to cut down on code clutter


def paramRecord(open_csv, starname, atmparams):
    open_csv.write(str(starname)+','+str(atmparams[0])+','+str(atmparams[1])+','+str(atmparams[2])+',' +
                   str(atmparams[3])+','+str(atmparams[4])+','+str(atmparams[5])+'\n')

# OKAY 05/25/16
# opens and parses a parameter file to postprocess


def paramParse(parameter_file):
    paramfile = open(parameter_file, 'r')
    paramlines = paramfile.readlines()
    smo_resolu = 0.0
    rottuple = (0.0, 0.0)
    rotfudge = 0.0
    radtuple = (0.0, 0.0)
    num_outputs = 0
    package = 'NO'
    verbose = 'NO'
    resolution = 0.00
    waverange = (0.0, 0.0)
    outspacing = 0.0
    for entry in paramlines:
        if 'NUMBER_OUTPUTS' in entry:
            num_outputs = entry.split(':')[1]
            num_outputs = int(num_outputs)
        if 'SMO_RESOLU' in entry:
            smo_resolu = entry.split(':')[1]
            smo_resolu = float(smo_resolu.strip())
        if 'ROT_VELOCI' in entry:
            rottuple = entry.split(':')[1]
            rottuple = rottuple.split(',')
            rottuple = (float(rottuple[0].rstrip()), float(rottuple[1].rstrip()))
        if 'ROT_FUDGE' in entry:
            rotfudge = entry.split(':')[1]
            rotfudge = float(rotfudge.rstrip())
        if 'RAD_VELOCI' in entry:
            radtuple = entry.split(':')[1]
            radtuple = radtuple.split(',')
            radtuple = (float(radtuple[0].rstrip()), float(radtuple[1].rstrip()))
        if 'BINARY_OUTPUT' in entry:
            package = entry.split(':')[1]
            if package.strip() == 'YES':
                package = 'YES'
            else:
                package = 'NO'
        if 'NAT_RESOLU' in entry:
            resolution = entry.split(':')[1]
            resolution = float(resolution.rstrip())
        if 'VERBOSE_OUTPUT' in entry:
            verbose = entry.split(':')[1]
            if verbose.strip() == 'YES':
                verbose = 'YES'
            else:
                verbose = 'NO'
        if 'WAVE_RANGE' in entry:
            waverange = entry.split(':')[1]
            waverange = waverange.split(',')
            waverange = (float(waverange[0].rstrip()), float(waverange[1].rstrip()))
        if 'PIXEL_SCALE' in entry:
            outspacing = float(entry.split(':')[1])
    if smo_resolu == 0.0:
        smo_resolu = 0.01
        print('warning, no output resolution specified, default to 0.01 A/px')
    if rottuple == (0.0, 0.0):
        rottuple = (0.0, 5.0)
        print('warning, no rotational velocity specified, default to [0.0, 5.0] km/s')
    if rotfudge == 0.0:
        rotfudge = 0.0
        print('warning, no rotational velocity fudge factor specified, default to 0.0 km/s')
    if num_outputs == 0:
        num_outputs = 1000
        print('warning, specify number of output spectra, default to 1000')
    if resolution == 0.00:
        resolution = 0.01
        print('warning, input resolution not specified, defaulting to 0.01')
    if waverange == (0.0, 0.0):
        print('warning, no output wavelength range specified, default to input range')
    if outspacing == 0.0:
        print('warning, no output pixel scale specified, default to 0.2038 A/px')
    paramfile.close()
    return rottuple, smo_resolu, radtuple, resolution, waverange, outspacing, rotfudge, package, num_outputs, verbose
# adds macroturbulent velocity to the spectrum

""" No longer needed, added to specgen v1.6+
def macTurb(starname, starname_output, spacing, velocity=3.8):
    subprocess.call(['macturb', starname, starname_output+'MT', str(spacing), str(velocity)])
"""

# OKAY 05/26/16
# adds in rotational broadening to a spectrum
# most of the things are hardcoded for now
# 04/27/17 added in a fudge factor to change the rotation by an offset

def rotVet(starname, starname_output, rotation_range, spacing, fudge):
    rotation = random.uniform(rotation_range[0], rotation_range[1])
    rotation = rotation + fudge
    subprocess.call(['avsini', starname, starname_output+'R', str(rotation), '0.6', str(spacing)])
    rotation = round(rotation - fudge, 1)
    return rotation


# OKAY 05/26/16
# smooths spectrum to target resolution


def smoothSpec(starname, inspacing, outresolution, outspacing, final_output_name='.vsn.spec'):
    output_name = starname.rstrip('.spec')
    subprocess.call(['smooth2', starname+'RS', output_name+final_output_name, str(inspacing), str(outresolution),
                     str(outspacing)])
    subprocess.call(['rm', starname+'RS'])


def radVel(starname, rv_range):
    """Adds radial velocity shift to a star based on a tuple of ranges (in km/s)
    Does this by shifting the wavelength grid and remapping the original flux values to that shifted grid,
    Then interpolating those flux values back into the original wavelength grid"""
    starname = str(starname)
    spectrum = pd.read_csv(starname+'R', header=None, delim_whitespace=True)
    initial_flux = spectrum[1].values  # fetch the unshifted flux values
    wavelength = spectrum[0].values  # fetch the initial wavelength grid
    radvel = random.uniform(rv_range[0], rv_range[1])  # draw the radvel
    shifted_wave = wavelength + (radvel*wavelength)/3.0e5  # calculate the shifted wavelength grid
    shifted_flux = ip.interp1d(shifted_wave, initial_flux, kind='linear', bounds_error=False, fill_value=1.0)
    spectrum[1] = shifted_flux(wavelength)  # interpolate shifted spectrum onto initial wavelength grid
    spectrum.to_csv(starname+'RS', sep=' ', header=False, index=False)
    subprocess.call(['rm', starname+'R'])
    return radvel

"""
# OKAY 05/26/16
# changes the signal to noise ratio from infinite to some value, using 1.0 as the signal level
# using gaussian noise for now, can change later. UPDATE 01/26/17: Now no longer adds signal-to-noise,
# just picks an SN value


def signalNoise(starname, snr_range, final_output_name='.vsn.spec'):
    starname = str(starname)
    spectrum = pd.read_csv(starname+'RS', header=None, delim_whitespace=True)
    spectrum_length = len(spectrum[0].values)
    snr = random.uniform(snr_range[0], snr_range[1])
    snr = round(snr, 1)
    stddev = 1.0/snr
    noise_vector = pd.DataFrame(np.random.normal(0.0, stddev, (spectrum_length, 1)), index=range(spectrum_length),
                                columns=['noise'])
    spectrum[1] = spectrum[1]+noise_vector['noise']
    output_name = starname.rstrip('.spec')
    spectrum.to_csv(output_name+final_output_name, sep=' ', header=False, index=False)
    subprocess.call(['rm', starname+'RS'])
    return snr
"""


def specTrim(starname, waverange, final_output_name='.vsn.spec'):
    """Trims an output file (*.vsn.spec) to specified wavelength range
    Honestly not necessary and just churns the disk right now, but
    helps with code modularity/clarity
    """
    starname = str(starname)
    output_name = starname.rstrip('.spec')
    spectrum = np.genfromtxt(output_name+final_output_name)
    spectrum = spectrum[(waverange[0] <= spectrum[:, 0]) & (spectrum[:, 0] <= waverange[1]), :]
    np.savetxt(output_name+final_output_name, spectrum)
# OKAY 05/26/16
# Implement a worker to do multiprocessing
# takes as input the tuple returned by paramParse (minus the numstars designation), returns the atm params to put in
# paramRecord


def worker(starparams, inputparams):
    starname = str(int(starparams[1]))+'.spec'
    starname_output = str(starparams[1]+(starparams[0]))+'.spec'
    spparams = starparams[2:]
    waverange = inputparams[4]
    # post-processing: rotational broadening, smoothing, macroturbulent
    rotvet = rotVet(starname, starname_output, inputparams[0], inputparams[3], inputparams[6])
    rv = radVel(starname_output, inputparams[2])
    smoothSpec(starname_output, inputparams[3], inputparams[1], inputparams[5])
    specTrim(starname_output, waverange)
    # returns the parameters that were used
    atmparams = (str(spparams[0]), str(spparams[1]), str(spparams[2]), str(spparams[3]), str(rotvet), str(rv))
    return starname_output, atmparams

# reads in an input file, and associated atmospheric parameters, and records output into a binary form.
# output file will be of the form described in the code preamble


def binaryPack(input_file, output_file, atm_parameters, first='yes'):
    input_spec = np.genfromtxt(input_file)
    if first == 'yes':
        wavelength = input_spec[:, 0]
        # dumps #px, #params, wavelengths to beginning of file
        pwavelength = np.zeros(len(wavelength) + 2)
        pwavelength[0] = float(len(wavelength))
        pwavelength[1] = float(len(atm_parameters)-2)
        pwavelength[2:] = wavelength
        w_array = array('d', pwavelength)
        # now params, fluxes for first star
        a_array = array('d', atm_parameters)
        flux = input_spec[:, 1]
        f_array = array('d', flux)
        w_array.tofile(output_file)
        a_array.tofile(output_file)
        f_array.tofile(output_file)
        # subprocess.call(['rm', input_file])
    else:
        # now just dumps params, fluxes for stars 2 to N
        flux = input_spec[:, 1]
        f_array = array('d', flux)
        a_array = array('d', atm_parameters)
        a_array.tofile(output_file)
        f_array.tofile(output_file)
        # subprocess.call(['rm', input_file])

# OKAY 05/26/16
# finally, basic control flow


def Main(parameterfile):
    # initialize time for factoring the code
    start_time = time.time()
    # open and extract parameters from parameter file
    allparams = paramParse(parameterfile)
    # read in parameters used to generate the spectra
    spectra_params = pd.read_csv('atmparams.out')
    print('constructing output file list')
    # shuffle the read-in parameters
    # separate out the postprocessing parameters
    postproc_params = allparams[:7]
    # generate list of files to be postprocessed
    num_highres = len(spectra_params.index)
    num_outputs = allparams[8]
    partial_set = num_outputs % num_highres
    num_repetitions = int(math.floor(num_outputs / num_highres))
    total_params = pd.DataFrame()
    if partial_set != 0:
        part_rep = pd.DataFrame(np.full((partial_set, 1), float(0.0000)))
        part_data = spectra_params[:partial_set]
        part_params = pd.concat([part_rep, part_data], axis=1)
        total_params = total_params.append(part_params, ignore_index=True)
    if num_repetitions > 0:
        for i in range(num_repetitions):
            temp_params = spectra_params
            temp_repetition = pd.DataFrame(np.full((num_highres, 1), float(i+1)*0.0001))
            temp_concat = pd.concat([temp_repetition, temp_params], axis=1)
            total_params = total_params.append(temp_concat)
    shuffled_params = (total_params.as_matrix())
    np.random.shuffle(shuffled_params)
    shuffled_params = tuple(map(tuple, shuffled_params))
    print('done; starting postprocessing')
    # opens output file
    outputfile = open('atmparams.postprocess.out', 'w')
    outputfile.write('starname,temp,grav,met,vt,rot,rv\n')
    # defines a partial function because pool.map only takes a single list of iterables
    # the partial function basically freezes the inputparams variable in worker to be modelparams
    # (for any instance of partial_worker, inputparams is always modelparams - turn into a constant
    partial_worker = partial(worker, inputparams=postproc_params)
    # initializes the pool of processes
    pool = multiprocessing.Pool(7)
    # maps partial_worker and list of stars to the pool, stores used parameters in a list
    outputs = pool.map(partial_worker, shuffled_params)
    # end the list of functions to go to pool
    pool.close()
    # wait for all processes to return
    pool.join()
    # record the used parameters
    print('done; recording output')
    for entry in outputs:
        starname = str(entry[0])
        starname = starname.rstrip('.spec')
        atmparams = entry[1]
        paramRecord(outputfile, starname, atmparams)
    outputfile.close()
    # total time to execute
    elapsed_time = time.time() - start_time
    print("Total Execute Time: " + str(elapsed_time) + "s")
    # binary output if desired (it's on by default)
    if allparams[7] == 'YES':
        print('binary output construction')
        # open binary file
        binary_output = open('allspec', 'wb')
        atmparams = np.genfromtxt('atmparams.postprocess.out', delimiter=',', skip_header=1)
        # write everything to binary file
        for i in range(len(atmparams[:, 0])):
            starname = str(atmparams[i, 0])
            if i == 0:
                binaryPack(starname+'.vsn.spec', binary_output, atmparams[i, :], first='yes')
            else:
                binaryPack(starname+'.vsn.spec', binary_output, atmparams[i, :], first='no')
        binary_output.close()
        # remove *.vsn.spec files if specified (by default this is on)
        if allparams[9] == 'NO':
            filelist = glob.glob('*.vsn.spec')
            print('binary output only; removing individual postprocessed spectra')
            for fl in filelist:
                subprocess.call(['rm', fl])

if __name__ == "__main__":
    Main('postproc.param')
