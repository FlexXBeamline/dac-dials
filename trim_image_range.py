"""
Trim image range based on mean I/sigmaI cutoff.

Notes:
- Reflections are removed from the beginning & end of the rotation dataset (the .expt file is not modified)
"""
import logging
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

import iotbx.phil
from dxtbx.model import ExperimentList
from dxtbx import flumpy

from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

help_message = __doc__

phil_scope = iotbx.phil.parse(
"""\
  threshold = 1.0
    .type = float
    .help = "Frames are trimmed from beginning/end if the mean I/sigma_I is less than this threshold"
  window = 10
    .type = int
    .help = "Moving window (number of images) for smoothing results, helpful especially for fine-sliced or sparse data"
    
  output {
     reflections = 'trimmed.refl'
        .type = str
        .help = "The trimmed reflections output filename"
     log = 'trim_image_range.log'
        .type = str
        .help = "Name of log file"
    }
""",
    process_includes=True,
)

logger = logging.getLogger('dials.command_line.trim_image_range')

def ffill(arr):
    mask = np.isnan(arr)
    idx = np.where(~mask,np.arange(mask.size),0)
    np.maximum.accumulate(idx,axis=0, out=idx)
    mask = np.logical_and(mask, idx)
    out = arr.copy()
    out[mask] = arr[idx[mask]]
    return out

def bfill(arr):
    return ffill(arr[::-1])[::-1]

def accum_sliding_window(v, winsize):
    vp = np.pad(v,((winsize-1)//2,winsize//2), constant_values=0)
    vw = sliding_window_view(vp, window_shape=winsize)
    va = vw.sum(axis=1)
    return va

def image_range_by_i_sig_i_cutoff(refl, threshold = 1.0, winsize = 10, minpoints = 10):
    # refl -- integrated reflections. use refl.select(refl.get_flags(refl.flags.integrated))

    # calculate i_sig_i, z for each reflection
    sel = refl['intensity.sum.variance'] > 0
    a = refl['intensity.sum.value'].select(sel).as_numpy_array()
    b = refl['intensity.sum.variance'].select(sel).as_numpy_array()
    i_sig_i = a/np.sqrt(b)
    z = refl['xyzcal.px'].select(sel).parts()[2].as_numpy_array()

    # bin by image number
    imgnum = np.floor(z).astype(int)
    s = np.bincount(imgnum, weights = i_sig_i)
    d = np.bincount(imgnum)
    x = range(0, 1 + np.max(imgnum))

    # smooth by sliding window

    sa = accum_sliding_window(s, winsize)
    da = accum_sliding_window(d, winsize)

    # temporarily suppress divide by zero warning
    np.seterr(divide='ignore', invalid='ignore') 
    y = bfill(ffill(sa/da))
    np.seterr(divide='warn', invalid='warn') # Reset to default
    ind = np.where(np.logical_and(y > threshold, da > minpoints))[0]
    
    return ind[0], ind[-1]

def trim_reflections_to_image_range(refl, zmin, zmax):
    z = refl['xyzcal.px'].parts()[2].as_numpy_array()
    isinrange = np.logical_and(z >= zmin, z < (zmax + 1))
    trimmed = refl.select(flumpy.from_numpy(isinrange))
    return trimmed

def run(args=None):
    
    usage = 'dials.python trim_image_range.py integrated.expt masked.refl [options]'
    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        epilog=help_message,
        check_format=False,
    )

    params, options, args = parser.parse_args(
        args, show_diff_phil=False, return_unhandled=True,
    )
    
    log.config(verbosity=options.verbose, logfile=params.output.log)
    
    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)
    
    refl, expt = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments,
    )
    refl = refl[0] # why is this a list?
    expt = expt[0]

    # run the algorithm
    integrated = refl.select(refl.get_flags(refl.flags.integrated))
    zmin, zmax = image_range_by_i_sig_i_cutoff(integrated, threshold = params.threshold, winsize=params.window)
    trimmed = trim_reflections_to_image_range(refl, zmin, zmax)

    # summary to log
    logger.info('-'*80)
    logger.info(f'trimmed image range to: {zmin}, {zmax}')
    logger.info(f'number of reflections: {len(refl)} (input) --> {len(trimmed)} (trimmed)')
    logger.info(f'Saving remaining reflections to {params.output.reflections}')
    logger.info('-'*80)
    
    # save the data
    trimmed.as_file(params.output.reflections)
    
if __name__ == "__main__":
    run()