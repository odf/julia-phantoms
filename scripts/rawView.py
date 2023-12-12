#!/usr/bin/env python

import sys
import numpy as np
from scipy import misc

filename = sys.argv[1] 
data = np.fromfile(filename,dtype=np.float32)
size = int(np.sqrt(len(data)))
image = data.reshape((size,size))
misc.toimage(image).show()

exit()
