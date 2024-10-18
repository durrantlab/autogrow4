# !/usr/bin/env python

"""This is the executable file for Autogrow 4.0.3. This script should come
first. It should obtain and verify all the parameters work. This than should
pass these parameters variables to the main execution function titled
AutogrowMainExecute.py found in MainFunctions

If you use AutoGrow 4.0.3 in your research, please cite the following reference:
Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm
for de novo drug design and lead optimization. J Cheminform 12, 25 (2020).
[doi: 10.1186/s13321-020-00429-4]
"""


import __future__
from autogrow.config.argparser import get_user_params  # TODO: It is strange that this is in argparser. Good to clean up user_params.py and move it there.
import autogrow.main as AutogrowMainExecute


################
# Run AutoGrow #
################

if __name__ == "__main__":
    params = get_user_params()
    AutogrowMainExecute.main(params)
