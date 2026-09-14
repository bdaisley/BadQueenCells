#!/usr/bin/env python3

import sys
import hydraMPP
print("Version:", hydraMPP.__version__)
print("Location:", hydraMPP.__file__)
print("First sys.path entries:")
print("\n".join(sys.path[:5]))