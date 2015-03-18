#!/usr/bin/env python

import numpy as np
from scipy.signal import lfilter, tf2zpk


b = np.array([1.0770389918072, 2.03115724017531, 0.961079950480885])
a = np.array([1.0000000000000, 1.98606828810359, 0.992414110487639])

