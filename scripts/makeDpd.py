#!/usr/bin/env python

import datetime
import glob
import optparse
import logging
import os
import re
import sys
from collections import namedtuple

import numpy as np

import mango


logger, rootLogger = mango.mpi.getLoggers(__name__)


def readPixelSize(path):
    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('pixelSize:'):
                return float(line.split()[1])


ProjectionGeometries = namedtuple(
    'ProjectionGeometries',
    ['sourcePos', 'detectorPos', 'wStep', 'hStep'])


def readProjectionGeometries(path):
    data = np.loadtxt(path, delimiter=',')

    sourcePos = (data[:, 2::-1]).astype(np.float32)
    detectorPos = (data[:, 5:2:-1]).astype(np.float32)
    wStep = (data[:, 8:5:-1]).astype(np.float32)
    hStep = (data[:, 11:8:-1]).astype(np.float32)

    return ProjectionGeometries(sourcePos, detectorPos, wStep, hStep)


def writeProjectionGeometries(filename, pg):
    with open(filename, "w") as fp:
        fp.write("# source_x,   source_y,   source_z,")
        fp.write("     detector_x,   detector_y,   detector_z,")
        fp.write("      wstep_x,    wstep_y,    wstep_z,")
        fp.write("      hstep_x,    hstep_y,    hstep_z\n")

        xyz = lambda a: tuple(a[::-1])

        for i in range(len(pg.sourcePos)):
            fp.write("%10.7f, %10.7f, %10.7f,   " % xyz(pg.sourcePos[i]))
            fp.write("%12.7f, %12.7f, %12.7f,   " % xyz(pg.detectorPos[i]))
            fp.write("%10.7f, %10.7f, %10.7f,   " % xyz(pg.wStep[i]))
            fp.write("%10.7f, %10.7f, %10.7f\n" % xyz(pg.hStep[i]))


def fmdListFromProjGeoms(projGeoms, active):
    dist = lambda vs: np.sqrt(vs[:,1]**2 + vs[:,2]**2)

    t = np.arctan2(-projGeoms.sourcePos[:, 1], -projGeoms.sourcePos[:, 2])
    sampleTheta = (t / np.pi * 180.0) % 360.0

    sampleZ = projGeoms.sourcePos[:, 0]
    sampleY = dist(projGeoms.sourcePos)
    cameraY = dist(projGeoms.detectorPos - projGeoms.sourcePos)

    fmdList = []

    for i in range(len(projGeoms.sourcePos)):
        if i in active:
            fmd = mango.recon.FrameMetaData()
            fmd.index = i
            fmd.sampleTheta = sampleTheta[i]
            fmd.sampleZPos = sampleZ[i]
            fmd.sampleYPos = sampleY[i]
            fmd.cameraYPos = cameraY[i]
            fmd.temperature = 19.0
            currtime = datetime.datetime.fromtimestamp(i)
            fmd.time = currtime.strftime("%Y%m%d_%H%M%S")
            fmdList.append(fmd)

    return fmdList


def rawData(path):
    return np.fromfile(path, dtype=np.float32)


def readProjections(projectionFiles, startidx=0, numproj=-1):
    if numproj <= 0:
        numproj = len(projectionFiles)

    data = rawData(projectionFiles[0])
    size = int(np.sqrt(len(data)))

    projDdsShape = (numproj, size, size)
    projDdsOrigin = (startidx, -(size // 2), -(size // 2))

    projDds = mango.zeros(
        shape=projDdsShape,
        mpidims=(0, 1, 1),
        origin=projDdsOrigin,
        halo=(0, 0, 0),
        dtype=np.float32,
        mtype="tomo_float")

    begin = projDds.subd.origin[0]
    end = begin + projDds.subd.shape[0]

    for i, path in enumerate(projectionFiles):
        if begin <= i < end:
            img = rawData(path).reshape(size, size)
            projDds.subd.asarray()[i - begin] = img

    return projDds


def run():
    parser = optparse.OptionParser("usage: %prog [OPTIONS] INDIR [TEMPLATE]")

    (options, args) = parser.parse_args()

    indir = re.sub(r'/$', '', args[0])
    name = os.path.basename(indir)

    if len(args) > 1:
        dpdPath = re.sub(r'/$', '', args[1])
        dpd = mango.recon.io.readDpd(dpdPath)
        pixelSz = dpd.md.getVoxelSize("mm")[1]
        geom = dpd.getGeometry()
    else:
        pixelSz = readPixelSize(os.path.join(indir, 'info.dat'))
        geom = mango.recon.ConeBeamGeometry()
        geom.readFromFile(os.path.join(indir, 'geom.dat'))

    projGeoms = readProjectionGeometries(os.path.join(indir, 'proj.dat'))
    writeProjectionGeometries("geoms_%s.dat" % name, projGeoms)

    globPattern = os.path.join(indir, 'proj-[0-9][0-9]*.raw')
    projectionFiles = sorted(glob.glob(globPattern))

    regexPattern = os.path.join(indir, 'proj-([0-9][0-9]*).raw')
    active = [
        int(re.sub(regexPattern, r'\1', path)) - 1
        for path in projectionFiles
    ]
    fmd = fmdListFromProjGeoms(projGeoms, active)

    numproj = len(projectionFiles)
    data = rawData(projectionFiles[0])
    size = int(np.sqrt(len(data)))

    shape = (numproj, size, size)
    origin = (0, -(size // 2), -(size // 2))

    outDpdPath = mango.recon.io.zeros(
        "projf32_%s.nc" % name, p_fmd=fmd,
        shape=shape, origin=origin,
        voxel_size=(1, pixelSz, pixelSz), voxel_unit='mm'
    )

    # Mango crashes when writing Dpd subdomains with n < 2 projections.
    # So we make all subdomains approximately equal in size, 500 <= n <= 1000.

    nrChunks = int(np.ceil(numproj / 1000.0))
    factor = (numproj + 0.001) / nrChunks
    offset = 0

    for i in range(nrChunks):
        start = int(i * factor)
        end = int((i + 1) * factor)

        chunkDds = readProjections(projectionFiles, start, end - start)
        chunkDds.md.setVoxelSize((1, pixelSz, pixelSz), "mm")

        rootLogger.info(
            "writing %s projections from %s" % (end - start, start)
        )

        chunkDpd = mango.recon.dpd(geom, chunkDds, fmd[start : end])
        mango.recon.io.writeDpdSubDomain(outDpdPath, chunkDpd)


if __name__ == "__main__":
    mango.setLoggingVerbosityLevel("high")
    mango.mpi.initialiseLoggers([__name__, "mango"], logLevel=logging.DEBUG)

    run()
