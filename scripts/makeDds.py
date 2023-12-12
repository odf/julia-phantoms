#!/usr/bin/env python

import glob
import optparse
import logging
import os
import re
import sys

import numpy as np

import mango


logger, rootLogger = mango.mpi.getLoggers(__name__)


def readVoxelSize(path):
    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('voxelSize:'):
                return float(line.split()[1])


def readVolumeSize(path):
    with open(path, 'r') as f:
        for line in f.readlines():
            if line.startswith('sizeXYZ:'):
                return map(float, line.split()[1:])


def rawData(path):
    return np.fromfile(path, dtype=np.float32)


def readSlices(sliceFiles, voxelSize, volumeSizeXYZ):
    numSlices = len(sliceFiles)
    xSize, ySize, _ = volumeSizeXYZ

    volDdsShape = (numSlices, ySize, xSize)
    volDdsOrigin = (-(numSlices // 2), -(ySize // 2), -(xSize // 2))

    volDds = mango.zeros(
        shape=volDdsShape,
        mpidims=(0, 1, 1),
        origin=volDdsOrigin,
        halo=(0, 0, 0),
        dtype=np.float32,
        mtype="tomo_float")

    volDds.md.setVoxelSize((voxelSize, voxelSize, voxelSize), "mm")

    begin = volDds.subd.origin[0] - volDds.origin[0]
    end = begin + volDds.subd.shape[0]

    for i, path in enumerate(sliceFiles):
        if begin <= i < end:
            img = rawData(path).reshape(ySize, xSize)
            volDds.subd.asarray()[i - begin] = img

    return volDds


def run():
    parser = optparse.OptionParser("usage: %prog [OPTIONS] INDIR")

    (options, args) = parser.parse_args()

    indir = re.sub(r'/$', '', args[0])
    name = os.path.basename(indir)

    voxelSz = readVoxelSize(os.path.join(indir, 'info.dat'))
    volumeSz = readVolumeSize(os.path.join(indir, 'info.dat'))

    pattern = os.path.join(indir, 'volz-[0-9][0-9]*.raw')
    sliceFiles = sorted(glob.glob(pattern))
    volDds = readSlices(sliceFiles, voxelSz, volumeSz)

    mango.io.writeDds("tomo_float_%s.nc" % name, volDds)


if __name__ == "__main__":
    mango.setLoggingVerbosityLevel("high")
    mango.mpi.initialiseLoggers([__name__, "mango"], logLevel=logging.DEBUG)

    run()
