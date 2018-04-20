"""
v1.4: this strips out all the post-processing and simply outputs infinite SN spectra at the specified resolution
(suggest 0.01). Each star is output as one file. A library of stars created using this script can then be
post-processed more quickly, enabling larger training set creation. NOTE: with resolution 0.01, 400 A wavelength range,
each file is approximately 850 kb.

Inputs: parameter file ('generate.param'), must have the ./atm kurucz directory in working directory,
linelist ('luke.extend.lst')

Outputs: *.spec files with columns wavelength, continuum-normalized flux, and summary of atm params ('atmparams.out')
"""


# reads in a parameter file that then gets parsed, and batch files suitable for generating Kurucz atm models and
# then feeding them into SPECTRUM with desired parameters are read out
# example parameter file is included in this directory
# works by randomly generating stars based on ranges of parameter inputs

import subprocess
import random
import multiprocessing
from functools import partial
import time
import os.path
import numpy as np
from array import array

# OKAY 05/24/16
# atmcreate takes in a starname, and lists/tuples with lower [0] and upper [1] bounds for atmospheric parameters
# generates an atmospheric model by making an external call, assuming the usual "atm" folder is in the working directory
# returns a tuple containing the the atmosphere filename and atmospheric parameters
# set up right now to do select gravities and vturbs with resolution 0.01, temperatures with resolution 1K
# and metallicities with resolution e[fe/h]=0.01 (does this by rounding, ignore the small edge bias this introduces)


def atmCreate(starname, temprange, gravrange, metalrange, vturbrange):
    atmname = 'tempatm.' + str(starname)
    temp = random.uniform(temprange[0], temprange[1])
    temp = round(temp, 0)
    grav = random.uniform(gravrange[0], gravrange[1])
    grav = round(grav, 2)
    metal = random.uniform(metalrange[0], metalrange[1])
    metal = round(metal, 3)
    vturb = random.uniform(vturbrange[0], vturbrange[1])
    vturb = round(vturb, 2)
    if metal < 0.0:
        mflag = '-m'
        # kurucz models use -m to denote negative metallicities so trim out the negatives when making the kurucz command
        kmetal = abs(metal)
    else:
        mflag = '-p'
        kmetal = metal
    subprocess.call(["./atm/mspawn", '-w'+atmname, '-t'+str(temp), '-g'+str(grav), mflag+str(kmetal), '-v'+str(vturb)])
    atmparams = (temp, grav, metal, vturb)
    return atmname, atmparams

# OKAY 05/24/16
# gets the atmospheric models into the correct format for SPECTRUM (make sure input is tempatm.X)
# go ahead and pass it the relevant atmospheric parameters so I don't have to dig them out of the atm file


def atmMunge(atmname, atmparams):
    readatm = open(atmname, 'r')
    atmlines = readatm.readlines()
    nlayers = int(atmlines[2].split()[1])
    atmosphere = atmlines[3:3+nlayers]
    starname = str(atmname.split('.')[1])
    writeatm = open('specatm.'+starname, 'w')
    writeatm.write(str(atmparams[0])+' '+str(atmparams[1])+' '+str(atmparams[2])+' '+str(nlayers)+' \n')
    for entry in atmosphere:
        writeatm.write(entry)
    readatm.close()
    writeatm.close()
    return 'specatm.'+starname

# OKAY 05/24/16
# creates the SPECTRUM batch input files, since apparently it can't handle doing it in one line
# first two inputs are strings, vturb is float, wavelength is a tuple with (lower,upper) in angstroms,
# resolution is float in angstroms


def specCreate(specready_atmname, linelist, vturb, wavelength, resolution):
    readatm = open(specready_atmname, 'r')
    starname = str(specready_atmname.split('.')[1])
    specparams = open('specparams.'+starname, 'w')
    outputname = starname+'.speci'
    specparams.write(specready_atmname+'\n'+linelist+'\n'+outputname+'\n'+str(vturb)+'\n'+str(wavelength[0]) +
                     ','+str(wavelength[1])+'\n'+str(resolution)+'\n\n')
    specparams.close()
    readatm.close()
    return 'specparams.'+starname

# OKAY 05/24/16
# finally executes SPECTRUM in silent mode on a parameter file


def specExecute(specparams):
    # using shell=True here under the understanding that everything being passed is internally-generated
    # so security is not an issue
    subprocess.call(['spectrum n < '+specparams], shell=True)


#added macroturbulent velocity calculation 04/28/2017
#relation from Doyle et al. 2016

def macTurb(spready_atmname, params, resolution):
    starname = str(spready_atmname.split('.')[1])
    velocity = 3.21 + 2.33e-3*(params[0]-5777.0)+2.0e-6*(params[0]-5777.0)**2 - 2.00*(params[1]-4.44)
    subprocess.call(['macturb', starname+'.speci', starname+'.spec', str(resolution), str(velocity)])
    subprocess.call(['rm', starname+'.speci'])

# OKAY 05/24/16
# simply removes the temp files generated during the process
# if it significantly speeds up the program, I can use wildcards and run this on a bunch of temp files at once at
# various intervals in the code. will need to test. starname should be a string


def cleanUp(starname):
    starname = str(starname)
    subprocess.call(['rm', 'tempatm.'+starname])
    subprocess.call(['rm', 'specparams.'+starname])
    subprocess.call(['rm', 'specatm.'+starname])

# OKAY 05/25/16
# pass this function an open csv file,
# and it will write the starname and atmospheric parameters used in a comma separated way
# just doing this to cut down on code clutter


def paramRecord(open_csv, starname, atmparams):
    open_csv.write(str(starname)+','+str(atmparams[0])+','+str(atmparams[1])+','+str(atmparams[2])+',' +
                   str(atmparams[3])+'\n')

# OKAY 05/25/16
# opens and parses a parameter file to generate the atmospheres and spectra


def paramParse(parameter_file):
    paramfile = open(parameter_file, 'r')
    paramlines = paramfile.readlines()
    temptuple, gravtuple, mettuple, vttuple, wavetuple = ((0.0, 0.0) for i in range(5))
    resolution = 0.0
    numstars = 0
    linelist = 'None Specified'
    for entry in paramlines:
        if 'NUM_STARS' in entry:
            numstars = entry.split(':')[1]
            numstars = int(numstars.rstrip())
        if 'TEMP_RANGE' in entry:
            templist = entry.split(':')[1]
            temptuple = templist.split(',')
            temptuple = (float(temptuple[0].rstrip()), float(temptuple[1].rstrip()))
        if 'GRAV_RANGE' in entry:
            templist = entry.split(':')[1]
            gravtuple = templist.split(',')
            gravtuple = (float(gravtuple[0].rstrip()), float(gravtuple[1].rstrip()))
        if 'MET_RANGE' in entry:
            templist = entry.split(':')[1]
            mettuple = templist.split(',')
            mettuple = (float(mettuple[0].rstrip()), float(mettuple[1].rstrip()))
        if 'VT_RANGE' in entry:
            templist = entry.split(':')[1]
            vttuple = templist.split(',')
            vttuple = (float(vttuple[0].rstrip()), float(vttuple[1].rstrip()))
        if 'WAVE_RANGE' in entry:
            templist = entry.split(':')[1]
            wavetuple = templist.split(',')
            wavetuple = (float(wavetuple[0].rstrip()), float(wavetuple[1].rstrip()))
        if 'NAT_RESOLU' in entry:
            resolution = entry.split(':')[1]
            resolution = float(resolution.rstrip())
        if 'SLINELIST' in entry:
            linelist = entry.split(':')[1]
            linelist = linelist.strip()
    paramfile.close()
    return numstars, temptuple, gravtuple, mettuple, vttuple, wavetuple, resolution, linelist

# OKAY 05/26/16
# Implement a worker to do multiprocessing
# takes as input the tuple returned by paramParse (minus the numstars designation), returns the atm params to put in
# paramRecord


def worker(starname, inputparams):
    # randomizes the seed
    random.seed(int(time.time())+starname*17)
    np.random.seed(int(time.time())+starname*17)
    # creates atmosphere model using the input parameters and starname
    firststep = atmCreate(starname, inputparams[0], inputparams[1], inputparams[2], inputparams[3])
    # checks to make sure the kurucz model was created. if not, immediately returns the parameters and a flag
    if os.path.isfile(firststep[0]):
        # puts the kurucz output into a spectrum-readable format
        secondstep = atmMunge(firststep[0], firststep[1])
        # pulls out the atmospheric parameters from atmCreate
        atmparams = firststep[1]
        # creates the control file to be read into spectrum
        thirdstep = specCreate(secondstep, inputparams[6], atmparams[3], inputparams[4], inputparams[5])
        specExecute(thirdstep)
        #macroturbulent velocity calculation
        macTurb(spready_atmname=secondstep, params=firststep[1], resolution=inputparams[5])
        # cleans up garbage
        cleanUp(starname)
        # returns the parameters that were used
        atmparams = (str(atmparams[0]), str(atmparams[1]), str(atmparams[2]), str(atmparams[3]))
        return starname, atmparams
    else:
        defaultatm = firststep[1]
        default = (defaultatm[0], defaultatm[1], defaultatm[2], defaultatm[3])
        return 'atm_fault', default

# OKAY 05/26/16
# finally, basic control flow


def Main(parameterfile):
    # initialize time for factoring the code
    start_time = time.time()
    # open and extract parameters from parameter file
    allparams = paramParse(parameterfile)
    # number of stars to generate is given by the param file
    starlist = range(allparams[0])
    # separates out the desired parameter ranges
    modelparams = allparams[1:]
    # opens output file
    outputfile = open('atmparams.out', 'w')
    outputfile.write('starname,temp,grav,met,vt\n')
    # defines a partial function because pool.map only takes a single list of iterables
    # the partial function basically freezes the inputparams variable in worker to be modelparams
    # (for any instance of partial_worker, inputparams is always modelparams - turn into a constant
    partial_worker = partial(worker, inputparams=modelparams)
    # initializes the pool of processes
    pool = multiprocessing.Pool(7)
    # maps partial_worker and list of stars to the pool, stores used parameters in a list
    outputs = pool.map(partial_worker, starlist)
    # end the list of functions to go to pool
    pool.close()
    # wait for all processes to return
    pool.join()
    # record the used parameters
    for entry in outputs:
        starname = str(entry[0])
        atmparams = entry[1]
        paramRecord(outputfile, starname, atmparams)
    outputfile.close()
    # total time to execute
    elapsed_time = time.time() - start_time
    print("Total Execute Time: " + str(elapsed_time) + "s")

if __name__ == "__main__":
    Main('generate.param')
