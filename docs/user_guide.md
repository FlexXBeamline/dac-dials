# Guide to processing DAC data with `dac-dials`

This is intended as a general guide, but the examples assume data were collected at CHESS. The CHESS method is to load several crystals in the gasket at once, and collect ~50 degrees of rotation data from each one. With this wide rotation range, the beginning and ending frames are typically unusable because the beam is passing through the gasket material or otherwise blocked by the DAC frame.

## Organizing directories

First, create a directory for your project. If processing locally on a CHESS computer, this should be in the "aux" folder for your beamtime:

```
/nfs/chess/aux/cycles/<run-cycle>/id7b2/<investigator>/<date>/
```

where the names in brackets should be replaced with correct values. To keep the output files of DIALS straight, I recommend a directory structure like this:

```
.
├── dials_processing_notes.md
├── sample1
│   ├── ambient_crystal1
│   ├── ambient_crystal2
│   ├── ambient_crystal3
│   ├── 2kbar_crystal1
│   ├── 2kbar_crystal2
└── sample2
    ├── ambient_crystal1
    ├── ambient_crystal2
    ├── 2kbar_crystal1
    ├── 2kbar_crystal2
    ├── 3kbar_crystal1
    └── 3kbar_crystal2
```

I recommend keeping a log of the commands you run in a text file, as well as your notes about decisions made, so that you can easily reproduce processing later: in the example above, it's called `dials_processing_notes.md`.

Here, a "sample" refers to a particular kind of crystal (for instance, a space group or particular ligand soak), and the subfolders refer to single datasets collected on that crystal, possibly under different conditions such as pressure. For the purposes of this tutorial, I'll assume the above directory structure.

## Processing each dataset

For each sample, a typical DIALS workflow looks like this:

```bash
cd <dataset-directory>
dials.import <path-to-master-h5-file> [options]
dials.find_spots imported.expt
dials.index imported.expt strong.refl [options]
dials.refine indexed.expt indexed.refl
dials.integrate refined.expt refined.refl
dac.fit_gasket_mask integrated.expt integrated.refl
dac.trim_image_range integrated.expt masked.refl
dials.anvil_correction integrated.expt trimmed.refl [options]
```

Each step and various options are described below.

### 1. Importing (DIALS)

The command `dials.import` will extract geometry information stored in the `_master.h5` image and create a text (json) file, `imported.expt`. For the DAC experiments at 7b2 performed prior to November 2023, you'll have to override the detector distance, rotation angle, and spindle axis information because it was not written to the header. For data after November 2023, you only need to specify that the rotation axis is vertical, as follows:

```
cd sample1/ambient_crystal1
dials.import axis=0,1,0 <path-to-master-h5-file> 
```

On the CHESS network, the image path will be something like `/nfs/chess/raw/<run-cycle>/id7b2/<investigator>/<date>/<...>/<prefix>_master.h5`.

If import is successful, it will generate a file `imported.expt`. You can view information about the experiment like this:

```
dials.show imported.expt
```

This is helpful for making sure the geometry is correct. You can also inspect the images using:

```
dials.image_viewer imported.expt
```

### 2. Spotfinding, indexing, and integration (DIALS)

First, find all of the intense spots

```
dials.find_spots imported.expt
```

Then, index the reflections (find a lattice consistent with the data):

```
dials.index imported.expt strong.refl space_group=<sg>
```

Specifying the space group is optional, but for multi-crystal data, I've found it's useful because cell angle restraints can be used in geometry refinement. If indexing fails, this is probably a bad dataset and not worth processing further (but double-check the geometry -- see above).

Finally, refine the geometry using the indexed data:

```
dials.refine indexed.expt indexed.refl
```

The default parameters for integration usually work just fine.

```
dials.integrate refined.expt refined.refl d_min=<resolution>
```

The resolution cutoff `d_min` is optional -- if you know it already, or can estimate from the images, it's nice to apply it here because it can speed up the integration step substantially.

### 3. Masking the gasket shadow (`dac-dials`)

Next, we will fit a model for the gasket aperture.

```
dac.fit_gasket_mask.py integrated.expt integrated.refl
```

The script saves a file `masked.refl` where the integrated reflections that are occluded have been removed. 

Optional arguments can be specified if needed. To see the list of options, type:

```
dac.fit_gasket_mask -c -a1
```

The following output is produced:

```
Showing configuration parameters with:
  attributes_level = 1
  expert_level = 0

phi0 = 0
  .help = "value of rotation angle when beam passes through the gasket (center"
          "of scan)"
threshold = 4
  .help = "I/sigma_I threshold for reflections used to determine the mask"
          "boundary"
output {
  reflections = 'masked.refl'
    .help = "The masked reflections output filename"
  log = 'dac.fit_gasket_mask.log'
    .help = "Name of log file"
  gasket_mask = 'gasket_mask.json'
    .help = "The best gasket mask model"
}
input {
  experiments = None
    .help = "The experiment list file path"
  reflections = None
    .help = "The reflection table file path"
}
```

In particular, the program assumes that the diamond windows are normal to the beam when the rotation angle is zero. If not, the correct angle should be set using the option `phi0=<angle-in-degrees>`, and in the following step (`dials.anvil_correction`).

### 4. Trimming the image range

With the wide oscillation range typically used at CHESS, there will be frames at the beginning and end of data collection with little or not diffraction. These can be removed manually (for instance in `dials.import`), but for a more automatic approach you can run `dac.trim_image_range` next:

```bash
dac.trim_image_range integrated.expt masked.refl
```

A file `trimmed.refl` is saved. The default thresholds are designed to work in most cases, but if needed they can be adjusted using additional arguments, which can be listed as follows:

```
dac.trim_image_range -c -a1
```

The following output is produced:

```
Showing configuration parameters with:
  attributes_level = 1
  expert_level = 0

threshold = 1.0
  .help = "Frames are trimmed from beginning/end if the mean I/sigma_I is less"
          "than this threshold"
window = 10
  .help = "Moving window (number of images) for smoothing results, helpful"
          "especially for fine-sliced or sparse data"
output {
  reflections = 'trimmed.refl'
    .help = "The trimmed reflections output filename"
  log = 'dac.trim_image_range.log'
    .help = "Name of log file"
}
input {
  experiments = None
    .help = "The experiment list file path"
  reflections = None
    .help = "The reflection table file path"
}
```

### 5. Correcting diamond absorption

To run the anvil correction script included with DIALS, you have to specify the diamond thickness and the vector normal to the plane of the diamond window when the rotation angle is zero. At CHESS, the diamond plate on the "exit" side is 0.8 mm thick.

```
dials.anvil_correction integrated.expt trimmed.refl thickness=0.8 normal=0,0,1
```

The script saves a file `corrected.refl`. If the anvil correction is working properly, the R-factors should improve. I recommend processing `corrected.refl` and `trimmed.refl` separately to check whether the anvil correction was helpful or not.

## Merging multiple datasets

After you have performed the above steps for several isomorphous crystals, you can scale and merge them collectively. In this example, let's assume we're going to merge the sample1 data collected at 2kbar. The procedure is to go to the `sample1` directory, and then make a directory for processing and `cd` into it, and run `xia2.multiplex`. This program automatically merges multi-crystal data and outputs a lot of statistics and reports, which you should scrutinize. The merged and unmerged intensities will be output as an mtz file for further processing.

First, try merging without the anvil correction:

```bash
mkdir 2kbar_multiplex
cd 2kbar_multiplex
xia2.multiplex ../2kbar_*/trimmed.refl ../2kbar_*/integrated.expt
```

Then, try merging with the anvil correction:

```bash
mkdir 2kbar_multiplex_corrected
cd 2kbar_multiplex_corrected
xia2.multiplex ../2kbar_*/corrected.refl ../2kbar_*/integrated.expt
```

Compare the processing statistics for both the corrected and uncorrected data. A good place to start is the table called "Summary of merging statistics" in xia2.multiplex.log. You can also find more information in the various html reports.

## Tips / caveats

- A common issue with DAC data collection is low completeness, both because the rotation range is limited (50 degrees or so) and because the gasket tends to block high-resolution reflections at the beginning and end of the rotation range. Generally, to approach 100% completeness in all resolution shells, this usually means merging data from multiple crystals. Thus, a multi-crystal strategy is important to consider from the outset. However, for a multi-crystal strategy to succeed, you need to be certain that the crystals are truly isomorphous. In order to maximize the chance of having isomorphous crystals, we recommend loading many at once in the gasket. That way, the crystals are in the same environment and have experienced roughly the same pressure treatment. If this is not possible, then you may need to be especially careful in pressure measurement and timing of the experiments.

- The `dac-dials` scripts are experimental, and have only been tested on a handful of datasets. So always inspect the outputs critically (the DIALS html reports are a good place to start). If something looks wrong, please get in touch.