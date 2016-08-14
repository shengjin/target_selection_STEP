# myconfig.py:

import math
# NOTICE: The current version of this program requires:
#    crit_ra_radians = crit_dec_radians
#    AND
#    crit_ra_pm = crit_dec_pm


# Requirement of target star
crit_distance_pc = 15.0
crit_magnitude_v = 11.0

# Set the critical difference between the Proper Motions in Right Ascension
#     of the target and the reference stars. (mas/yr)
crit_pm_ratio = 0.01
crit_pm = False

# Set the critical difference between the Proper Motions in Declination
#     of the target and the reference stars. (mas/yr)
#crit_dec_pm = 10.0


##################
# para set by u4acess
#mag12, magnitude range of reference star
#ra0dc0, the epoch2000 of each target star
#wrawdc, the box-size for reference-selection
