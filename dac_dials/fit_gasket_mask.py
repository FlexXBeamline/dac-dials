"""
Fit a physical model for the region of the detector masked by the gasket during data collection
"""

import logging
import numpy as np

import iotbx.phil
from dxtbx.model import ExperimentList
from dxtbx import flumpy

from dials.array_family import flex
from dials.util import log
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

help_message = __doc__

phil_scope = iotbx.phil.parse(
    """\
  phi0 = 0
    .type = float
    .help = "value of rotation angle when beam passes through the gasket (center of scan)"
  outlierweight = 30
    .type = float
    .help = "relative weight for outliers vs. masked symmetric peaks in scoring function"
  maxoutliers = 50
    .type = int
    .help = "Maximum number of outliers to consider"
  i_sig_i_cutoff = 4
    .type = float
    .help = "I/sigma_I cutoff for reflections used to determine the mask boundary"
  maxpoints = 500
    .type = int
    .help = "Maximum number of points to consider (those kept are furthest from the mean aperture position)"
    
  output {
     reflections = 'masked.refl'
        .type = str
        .help = "The masked reflections output filename"
     log = 'dac.fit_gasket_mask.log'
        .type = str
        .help = "Name of log file"
    }
""",
    process_includes=True,
)

logger = logging.getLogger("dials.command_line.fit_gasket_mask")


def load_data(exptfile, reflfile):
    """read the experiment and indexed reflections"""
    expt = ExperimentList.from_file(exptfile)
    expt = expt[0]
    refl = flex.reflection_table.from_file(reflfile)
    refl = refl.select(refl.get_flags(refl.flags.indexed))
    return expt, refl


def calc_diffraction_vectors(expt, refl):
    panel = expt.detector[0]
    x, y, rot_angle = refl["xyzobs.mm.value"].parts()
    xyz = panel.get_lab_coord(flex.vec2_double(x, y))
    r = xyz / xyz.norms()
    r = flumpy.to_numpy(r)

    rotation_axis = expt.goniometer.get_rotation_axis()
    beam_direction = expt.beam.get_unit_s0()

    e1 = np.cross(beam_direction, rotation_axis)  # x
    e2 = np.cross(
        rotation_axis, e1
    )  # y, nominally along beam direction if rotation axis is perpendicular to beam
    e3 = rotation_axis  # z

    r1 = np.dot(e1, r.T)
    r2 = np.dot(e2, r.T)
    r3 = np.dot(e3, r.T)
    angle = flumpy.to_numpy(rot_angle) * (180 / np.pi)

    return angle, r1, r2, r3


def to_goniometer_frame(angle, r1, r2, r3, phi0=0):
    delta_phi = (np.pi / 180) * (angle - phi0)
    c, s = np.cos(delta_phi), np.sin(delta_phi)
    R11, R12 = c, -s
    R21, R22 = s, c

    # calculate r0 = R.T * r to rotate back into the frame of the spindle
    r1r = R11 * r1 + R21 * r2
    r2r = R12 * r1 + R22 * r2
    r3r = r3

    return r1r, r2r, r3r


def to_aperture_plane(r1, r2, r3):
    # x,y position of reflection in a plane parallel to the aperture plane, at a normal distance of "1 mm" from the crystal
    # or, equivalently, in the aperture plane in units where the normal distance to the crystal is "1"
    # b = r*(az/rz) - a
    # b = r*(1/r3) - [0,1,0]
    b1 = r1 / r2
    b3 = r3 / r2
    return b1, b3


class SmallestCircle:
    """find the smallest circle enclosing a set of points using the welzl algorithm"""

    def __init__(self, points):
        self.points = points

    def run(self, point_set=None, boundary_set=None):
        if point_set is None:
            point_set = tuple(range(self.points.shape[0]))
        else:
            point_set = tuple(point_set)
        if boundary_set is None:
            boundary_set = tuple()
        return self.welzl(point_set, boundary_set)

    def fit_circle_to_points(self, indices):
        def zero_points():
            return [0, 0], 0

        def one_point(p):
            return p, 0

        def two_points(p1, p2):
            center = p1 * 0.5 + p2 * 0.5
            r = p1 - center
            radius = np.sqrt(r[0] * r[0] + r[1] * r[1])
            return center, radius

        def three_points(p1, p2, p3):
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3
            A2 = x1 * x1 + y1 * y1
            B2 = x2 * x2 + y2 * y2
            C2 = x3 * x3 + y3 * y3
            Sx = 0.5 * np.linalg.det([[A2, B2, C2], [y1, y2, y3], [1, 1, 1]])
            Sy = 0.5 * np.linalg.det([[x1, x2, x3], [A2, B2, C2], [1, 1, 1]])
            a = np.linalg.det([[x1, x2, x3], [y1, y2, y3], [1, 1, 1]])
            b = np.linalg.det([[x1, x2, x3], [y1, y2, y3], [A2, B2, C2]])
            center = [Sx / a, Sy / a]
            radius = np.sqrt(b / a + (Sx * Sx + Sy * Sy) / (a * a))
            return center, radius

        points = self.points[indices, :]
        npoints = points.shape[0]
        if npoints == 0:
            return zero_points()
        elif npoints == 1:
            return one_point(points[0, :])
        elif npoints == 2:  # center at the mid-point
            return two_points(points[0, :], points[1, :])
        elif npoints == 3:  # center using triangle algorithm
            return three_points(points[0, :], points[1, :], points[2, :])

    def is_in_circle(self, center, radius, index):
        r = self.points[index, :] - center
        return r[0] * r[0] + r[1] * r[1] <= radius * radius

    def welzl(self, point_set, boundary_set):
        """the recursive function in welzl algorithm"""
        if len(point_set) == 0 or len(boundary_set) == 3:
            center, radius = self.fit_circle_to_points(boundary_set)
            return center, radius, boundary_set
        p = point_set[:1]
        center, radius, final_boundary = self.welzl(point_set[1:], boundary_set)
        if self.is_in_circle(center, radius, p[0]):
            return center, radius, final_boundary
        return self.welzl(point_set[1:], boundary_set + p)


def arg_furthest_points(x, y, Nmax, center=None):
    if Nmax > len(x):
        Nmax = len(x)
    if center is None:
        center = [np.mean(x), np.mean(y)]
    center_distance = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    return np.argpartition(center_distance, -Nmax)[-Nmax:]


def fit_circle(x, y, point_set=None, boundary_set=None, Nmax=1000):
    npoints = len(x)
    if point_set is not None:
        npoints = len(point_set)
    else:
        npoints = len(x)
    if npoints > Nmax:  # the max recursion depth is exceeded for too many points
        raise Exception("too many points, choose a subset to fit")
    sc = SmallestCircle(np.stack((x, y)).T)
    return sc.run(point_set=point_set, boundary_set=boundary_set)


def run(args=None):
    usage = "dials.python fit_gasket_mask.py integrated.expt integrated.refl [options]"
    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        epilog=help_message,
        check_format=False,
    )

    params, options, args = parser.parse_args(
        args,
        show_diff_phil=False,
        return_unhandled=True,
    )

    log.config(verbosity=options.verbose, logfile=params.output.log)

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    refl, expt = reflections_and_experiments_from_files(
        params.input.reflections,
        params.input.experiments,
    )
    refl = refl[0]  # why is this a list?
    expt = expt[0]

    # keep only indexed reflections if they exist
    # indexed = refl.select(refl.get_flags(refl.flags.indexed))

    integrated = refl.select(refl.get_flags(refl.flags.integrated))
    a = integrated["intensity.sum.value"].as_numpy_array()
    b = integrated["intensity.sum.variance"].as_numpy_array()
    # plt.plot(a/np.sqrt(b))
    isincl = a / np.sqrt(b) > params.i_sig_i_cutoff
    strong = integrated.select(flumpy.from_numpy(isincl))

    predicted = refl.select(refl.get_flags(refl.flags.predicted))

    # make sure we have indexed and predicted reflections present
    if strong.size() == 0 or predicted.size() == 0:
        logger.info(
            "Algorithm only works on integrated reflections. run dials.integrate first"
        )
        return
    logger.info(f"Keeping strong reflections only ({strong.size()} / {refl.size()})")
    refl = strong

    # convert to diffraction vectors
    # angle, r1, r2, r3 = calc_diffraction_vectors(expt, refl)
    b1, b3 = to_aperture_plane(
        *to_goniometer_frame(*calc_diffraction_vectors(expt, refl), phi0=params.phi0)
    )
    c1, c3 = to_aperture_plane(
        *to_goniometer_frame(
            *calc_diffraction_vectors(expt, predicted), phi0=params.phi0
        )
    )

    def remove_point(point_set, o):
        return tuple(set(point_set) - set((o,)))

    point_set = arg_furthest_points(b1, b3, params.maxpoints)
    boundary_set = None
    outlier_set = []
    results = []
    logger.info("Fitting models")
    for j in range(params.maxoutliers):
        center, radius, boundary_set = fit_circle(b1, b3, point_set=point_set)
        result = {
            "center": center,
            "radius": radius,
            "boundary_set": boundary_set,
            "outlier_set": tuple(outlier_set),
        }
        logger.info(
            f"center=[{center[0]:.5},{center[1]:.5}], radius={radius:7.5}, outliers={len(outlier_set)}"
        )
        results.append(result)
        bpindex = arg_furthest_points(
            b1[np.array(boundary_set)],
            b3[np.array(boundary_set)],
            1,
            center=[np.mean(b1), np.mean(b3)],
        )
        o = boundary_set[bpindex[0]]
        point_set = remove_point(point_set, o)
        boundary_set = remove_point(boundary_set, o)
        outlier_set.append(o)

    # compute scores
    scores = []
    for res in results:
        dx = c1 - res["center"][0]
        dy = c3 - res["center"][1]
        rsq = dx * dx + dy * dy
        num_outliers_calc = np.sum(rsq > (res["radius"] * res["radius"]))
        scores.append(
            num_outliers_calc - params.outlierweight * len(res["outlier_set"])
        )

    model = results[np.argmax(scores)]
    logger.info("-" * 80)
    logger.info("best-scoring model:")
    logger.info(
        f"center={model['center']}, radius={model['radius']}, outliers={len(model['outlier_set'])}"
    )
    logger.info("-" * 80)

    # apply the model

    # read the data again
    logger.info("reading reflection list again (is this necessary?)")
    refl, expt = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    refl = refl[0]  # why is this a list?
    expt = expt[0]

    logger.info("applying the model")
    angle, r1, r2, r3 = calc_diffraction_vectors(expt, refl)
    b1, b3 = to_aperture_plane(
        *to_goniometer_frame(angle, r1, r2, r3, phi0=params.phi0)
    )

    dx = b1 - model["center"][0]
    dy = b3 - model["center"][1]
    rsq = dx * dx + dy * dy
    isincl = rsq < 1.0001 * model["radius"] * model["radius"]

    masked = refl.select(flumpy.from_numpy(isincl))

    logger.info("-" * 80)
    logger.info(f"{masked.size()} reflections were kept (of {refl.size()})")
    logger.info(f"Saving remaining reflections to {params.output.reflections}")
    logger.info("-" * 80)

    # save the data
    masked.as_file(params.output.reflections)


if __name__ == "__main__":
    run()
