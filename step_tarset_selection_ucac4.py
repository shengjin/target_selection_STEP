#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

#Goto line 528 for the main program

"""
This script was written for the STEP space misssion to select target stars.
(c)  Sheng JIN, PMO, 2014

"""

try:
    import numpy as np
except:
    print ('ERROR: Numpy cannot be imported')
    print (' You need to inctall numpy ')
    quit()

try:
    import matplotlib.pyplot as plt
except:
    print ('WARNING:')
    print (' "matploblib.pylab" cannot be imported ')
    print (' To used the visualization functionality you need to install matplotlib')
    print (' Without matplotlib you can use this script for target selection ')
    print ('    but you will not be able to plot things or display images')


import os
import subprocess
import csv
import math
import time
import shutil
#import scipy
#import astropy
#import aplpy

try:
    from myconfig import *
except:
    print ('ERROR:')
    print ('The "myconfig.py" file is needed for the input parameters!')
    quit()


debug = True
#debug = False
#ddebug = True
ddebug = False


####################################
######## Class Defination
####################################

########## Class Hip2star 

class Hip2star():

    """
    Represents a single star in the Hipparcos new reduction data.
    """

    def __init__(self,
                 hipparcos_id,
                 solution_type_new,
                 solution_type_old,
                 number_of_components,
                 ra_radians,
                 dec_radians,
                 parallax_mas,
                 proper_motion_ra_mas_per_year,
                 proper_motion_dec_mas_per_year,
                 ra_error_mas,
                 dec_error_mas,
                 parallax_error_mas,
                 proper_motion_ra_error_mas_per_year,
                 proper_motion_dec_error_mas_per_year,
                 number_of_field_transits,
                 goodness_of_fit,
                 percentage_rejected_data,
                 cosmic_dispersion_added,
                 entry_in_supplemental_catalog,
                 magnitude,
                 magnitude_error,
                 magnitude_scatter,
                 variability_annex,
                 color_index,
                 color_index_error,
                 VI_color_index,
                 upper_triangular_weight_matrix):

        self.hipparcos_id = hipparcos_id
        self.solution_type_new = solution_type_new
        self.solution_type_old = solution_type_old
        self.number_of_components = number_of_components
        self.ra_radians = ra_radians
        self.dec_radians = dec_radians
        self.parallax_mas = parallax_mas
        self.proper_motion_ra_mas_per_year = proper_motion_ra_mas_per_year
        self.proper_motion_dec_mas_per_year = proper_motion_dec_mas_per_year
        self.ra_error_mas = ra_error_mas
        self.dec_error_mas = dec_error_mas
        self.parallax_error_mas = parallax_error_mas
        self.proper_motion_ra_error_mas_per_year = proper_motion_ra_error_mas_per_year
        self.proper_motion_dec_error_mas_per_year = proper_motion_dec_error_mas_per_year
        self.number_of_field_transits = number_of_field_transits
        self.goodness_of_fit = goodness_of_fit
        self.percentage_rejected_data = percentage_rejected_data
        self.cosmic_dispersion_added = cosmic_dispersion_added
        self.entry_in_supplemental_catalog = entry_in_supplemental_catalog
        self.magnitude = magnitude
        self.magnitude_error = magnitude_error
        self.magnitude_scatter = magnitude_scatter
        self.variability_annex = variability_annex
        self.color_index = color_index
        self.color_index_error = color_index_error
        self.VI_color_index = VI_color_index
        self.upper_triangular_weight_matrix = upper_triangular_weight_matrix

    def __str__(self):
        return( "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s" 
                 % (str(self.hipparcos_id), 
                    str(self.solution_type_new),
                    str(self.solution_type_old),
                    str(self.number_of_components), 
                    str(self.ra_radians),
                    str(self.dec_radians),
                    str(self.parallax_mas),
                    str(self.proper_motion_ra_mas_per_year),
                    str(self.proper_motion_dec_mas_per_year),
                    str(self.ra_error_mas),
                    str(self.dec_error_mas),
                    str(self.parallax_error_mas),
                    str(self.proper_motion_ra_error_mas_per_year),
                    str(self.proper_motion_dec_error_mas_per_year),
                    str(self.number_of_field_transits),
                    str(self.goodness_of_fit),
                    str(self.percentage_rejected_data),
                    str(self.cosmic_dispersion_added),
                    str(self.entry_in_supplemental_catalog),
                    str(self.magnitude),
                    str(self.magnitude_error),
                    str(self.magnitude_scatter),
                    str(self.variability_annex),
                    str(self.color_index),
                    str(self.color_index_error),
                    str(self.VI_color_index),
                    str(self.upper_triangular_weight_matrix)) )


############# Class HYGstar 

class HYGstar():

    """
    Represents a single star in the HYG (xyz) database.
    """

    def __init__(self,
                 hyg_id,
                 hipparcos_id,
                 hd_name,
                 hr_name,
                 gliese_name,
                 bayerflamsteed_name,
                 proper_name,
                 ra_epoch2000,
                 dec_epoch2000,
                 distance_pc,
                 pm_ra,
                 pm_dec,
                 rv,
                 magnitude_v,
                 absmagnitude_v,
                 spectrum,
                 color_index,
                 x_epoch2000,
                 y_epoch2000,
                 z_epoch2000,
                 v_x,
                 v_y,
                 v_z):

        self.hyg_id = hyg_id
        self.hipparcos_id = hipparcos_id
        self.hd_name = hd_name
        self.hr_name = hr_name
        self.gliese_name = gliese_name
        self.bayerflamsteed_name = bayerflamsteed_name
        self.proper_name = proper_name
        self.ra_epoch2000 = ra_epoch2000
        self.dec_epoch2000 = dec_epoch2000
        self.distance_pc = distance_pc
        self.pm_ra = pm_ra
        self.pm_dec = pm_dec
        self.rv = rv
        self.magnitude_v = magnitude_v
        self.absmagnitude_v = absmagnitude_v
        self.spectrum = spectrum
        self.color_index = color_index
        self.x_epoch2000 = x_epoch2000
        self.y_epoch2000 = y_epoch2000
        self.z_epoch2000 = z_epoch2000
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z

    def __str__(self):
        return( "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s" 
                 % (str(self.hyg_id),
                    str(self.hipparcos_id), 
                    str(self.hd_name),
                    str(self.hr_name),
                    str(self.gliese_name),
                    str(self.bayerflamsteed_name),
                    str(self.proper_name),
                    str(self.ra_epoch2000),
                    str(self.dec_epoch2000),
                    str(self.distance_pc),
                    str(self.pm_ra),
                    str(self.pm_dec),
                    str(self.rv),
                    str(self.magnitude_v), 
                    str(self.absmagnitude_v),
                    str(self.spectrum), 
                    str(self.color_index),
                    str(self.x_epoch2000),
                    str(self.y_epoch2000),
                    str(self.z_epoch2000),
                    str(self.v_x),
                    str(self.v_y),
                    str(self.v_z)) )



####################################
######## Function Defination
####################################


####################################
# Count the number of lines of a file
def count_lines(filename, numberline):
    """
    Count the number of lines of a file
    """
    numberline = sum(1 for line in open('filename')) 
    return numberline


####################################
# Check that if a file exists.
def check_file(filename, debug):
    """
    Check if the "filename" file exists.
    """
    if os.path.exists(filename):
        if debug:
            print ( 'The "%s" file exists.\n' % filename )
            print ()
    else:
        print ( 'The "%s" file does NOT exist! :o\n' % filename )
        quit()
    

####################################
# Read all the stars from the hygxyz.csv
def read_all_hyg_stars(filename):
    """
    Read all the information of all the HYG stars in "filename".
    """
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        star_list = [] 
        for row in reader:
            hyg_id = row[0]
            hipparcos_id = row[1]
            hd_name = row[2]
            hr_name = row[3]
            gliese_name = row[4]
            bayerflamsteed_name = row[5]
            proper_name = row[6]
            ra_epoch2000 = row[7]
            dec_epoch2000 = row[8]
            distance_pc = row[9]
            pm_ra = row[10]
            pm_dec = row[11]
            rv = row[12]
            magnitude_v = row[13]
            absmagnitude_v = row[14]
            spectrum = row[15]
            color_index = row[16]
            x_epoch2000 = row[17]
            y_epoch2000 = row[18]
            z_epoch2000 = row[19]
            v_x = row[20]
            v_y = row[21]
            v_z = row[22]

            star_once = HYGstar(hyg_id,
                        hipparcos_id,
                        hd_name,
                        hr_name,
                        gliese_name,
                        bayerflamsteed_name,
                        proper_name,
                        ra_epoch2000,
                        dec_epoch2000,
                        distance_pc,
                        pm_ra,
                        pm_dec,
                        rv,
                        magnitude_v,
                        absmagnitude_v,
                        spectrum,
                        color_index,
                        x_epoch2000,
                        y_epoch2000,
                        z_epoch2000,
                        v_x,
                        v_y,
                        v_z)
            star_list.append(star_once)
        return star_list


####################################
# Write all the hyg stars to the disk
def write_all_hyg_stars(star_list, filename):
    """
    Write all the information of all the HYG stars in "filename".
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['hyg_id',
                         'hipparcos_id',
                         'hd_name',
                         'hr_name',
                         'gliese_name',
                         'bayerflamsteed_name',
                         'proper_name',
                         'ra_epoch2000',
                         'dec_epoch2000',
                         'distance_pc',
                         'pm_ra',
                         'pm_dec',
                         'rv',
                         'magnitude_v',
                         'absmagnitude_v',
                         'spectrum',
                         'color_index',
                         'x_epoch2000',
                         'y_epoch2000',
                         'z_epoch2000',
                         'v_x',
                         'v_y',
                         'v_z'])

        for star in star_list:
            writer.writerow([star.hyg_id,
                             star.hipparcos_id,
                             star.hd_name,
                             star.hr_name,
                             star.gliese_name,
                             star.bayerflamsteed_name,
                             star.proper_name,
                             star.ra_epoch2000,
                             star.dec_epoch2000,
                             star.distance_pc,
                             star.pm_ra,
                             star.pm_dec,
                             star.rv,
                             star.magnitude_v,
                             star.absmagnitude_v,
                             star.spectrum,
                             star.color_index,
                             star.x_epoch2000, 
                             star.y_epoch2000,
                             star.z_epoch2000,
                             star.v_x, 
                             star.v_y,
                             star.v_z])


####################################
# Read all the stars from the hip2.dat
def read_all_hip2_stars(filename):
    """
    Read all the information of all the hip2 (Hipparcos 2007 A&A) stars in "filename".
    """
    with open(filename, "rb") as f:
        star_list = [] 
        row = f.read(277)
        while row:
            hipparcos_id = int(row[0:6])
            solution_type_new = int(row[7:10])
            solution_type_old = int(row[11])
            number_of_components = int(row[13])
            ra_radians = float(row[15:28])
            dec_radians = float(row[29:42])
            parallax_mas = float(row[43:50])
            proper_motion_ra_mas_per_year = float(row[51:59])
            proper_motion_dec_mas_per_year = float(row[60:68])
            ra_error_mas = float(row[69:75])
            dec_error_mas = float(row[76:82])
            parallax_error_mas = float(row[83:89])
            proper_motion_ra_error_mas_per_year = float(row[90:96])
            proper_motion_dec_error_mas_per_year = float(row[97:103])
            number_of_field_transits = int(row[104:107])
            goodness_of_fit = float(row[108:113])
            percentage_rejected_data = int(row[114:116])
            cosmic_dispersion_added = float(row[117:123])
            entry_in_supplemental_catalog = int(row[124:128])
            magnitude = float(row[129:136])
            magnitude_error = float(row[137:143])
            magnitude_scatter = float(row[144:149])
            variability_annex = int(row[150])
            color_index = float(row[152:158])
            color_index_error = float(row[159:164])
            VI_color_index = float(row[165:171])
            upper_triangular_weight_matrix = [[float(row[171:178]), float(row[178:185]), float(row[185:192]), float(row[192:199]), float(row[199:206])],
                                             [0.0, float(row[206:213]), float(row[213:220]), float(row[220:227]), float(row[227:234])],
                                             [0.0, 0.0, float(row[234:241]), float(row[241:248]), float(row[248:255])],
                                             [0.0, 0.0, 0.0, float(row[255:262]), float(row[262:269])],
                                             [0.0, 0.0, 0.0, 0.0, float(row[269:276])]]
    
            star_once = Hip2star(hipparcos_id,
                        solution_type_new,
                        solution_type_old,
                        number_of_components,
                        ra_radians,
                        dec_radians,
                        parallax_mas,
                        proper_motion_ra_mas_per_year,
                        proper_motion_dec_mas_per_year,
                        ra_error_mas,
                        dec_error_mas,
                        parallax_error_mas,
                        proper_motion_ra_error_mas_per_year,
                        proper_motion_dec_error_mas_per_year,
                        number_of_field_transits,
                        goodness_of_fit,
                        percentage_rejected_data,
                        cosmic_dispersion_added,
                        entry_in_supplemental_catalog,
                        magnitude,
                        magnitude_error,
                        magnitude_scatter,
                        variability_annex,
                        color_index,
                        color_index_error,
                        VI_color_index,
                        upper_triangular_weight_matrix)
            row = f.read(277)
            star_list.append(star_once)
        return star_list


####################################
# Write all the hyg stars to a csv file
def write_all_hip2_stars(star_list, filename):
    """
    Write all the information of all the hip2 (Hipparcos 2007 A&A) stars in "filename".
    WARNINT: Very strange behaver in the solution_type_new,old, and the variability_annex columns,
    WARNINT:      that is why I subtract 48 from these three values.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['hipparcos_id',
                         'solution_type_new',
                         'solution_type_old',
                         'number_of_components',
                         'ra_radians',
                         'dec_radians',
                         'parallax_mas',
                         'proper_motion_ra_mas_per_year',
                         'proper_motion_dec_mas_per_year',
                         'ra_error_mas',
                         'dec_error_mas',
                         'parallax_error_mas',
                         'proper_motion_ra_error_mas_per_year',
                         'proper_motion_dec_error_mas_per_year',
                         'number_of_field_transits',
                         'goodness_of_fit',
                         'percentage_rejected_data',
                         'cosmic_dispersion_added',
                         'entry_in_supplemental_catalog',
                         'magnitude',
                         'magnitude_error',
                         'magnitude_scatter',
                         'variability_annex',
                         'color_index',
                         'color_index_error',
                         'VI_color_index',
                         'upper_triangular_weight_matrix'])

        for star in star_list:
            writer.writerow([star.hipparcos_id,
                             star.solution_type_new,
                             star.solution_type_old-48,
                             star.number_of_components-48,
                             star.ra_radians,
                             star.dec_radians,
                             star.parallax_mas,
                             star.proper_motion_ra_mas_per_year,
                             star.proper_motion_dec_mas_per_year,
                             star.ra_error_mas,
                             star.dec_error_mas,
                             star.parallax_error_mas,
                             star.proper_motion_ra_error_mas_per_year,
                             star.proper_motion_dec_error_mas_per_year,
                             star.number_of_field_transits,
                             star.goodness_of_fit,
                             star.percentage_rejected_data,
                             star.cosmic_dispersion_added,
                             star.entry_in_supplemental_catalog,
                             star.magnitude,
                             star.magnitude_error,
                             star.magnitude_scatter,
                             star.variability_annex-48,
                             star.color_index,
                             star.color_index_error,
                             star.VI_color_index,
                             star.upper_triangular_weight_matrix])


####################################
######## The Main Program
####################################

if debug:
    print ()
    print ('#########################################')
    print ('  Target Selection for the STEP mission  ')
    print ('#########################################')
    print ()

### check if all the needed files exist
# the Hipparcos Catalog used to find refrence stars
#check_file('hip2.dat', debug)     
# the HYG Datebase used to find targets
check_file('hygxyz.csv', debug)


# Read all the stars from hygxyz.csv
if debug: print ('Reading the "hygxyz.cvs" file:')
hygstar_all = read_all_hyg_stars('hygxyz.csv')
hygstar_all_num = len(hygstar_all)
if debug: print ('    the HYG Database contains', hygstar_all_num, 'stars.\n')
if debug: print ()



# Select the hyg stars (F,G,K).
if debug: print ('Finding the F,G,K stars in the HYG database that is within', crit_distance_pc, 'pc and has a magnitude_v <', crit_magnitude_v, ':')
hygstar_selected = []
for i in range(hygstar_all_num):
#for i in range(100):
    distance_pc = float(hygstar_all[i].distance_pc)
    magnitude_v = float(hygstar_all[i].magnitude_v)
    if hygstar_all[i].spectrum:
        spectrum = hygstar_all[i].spectrum[0]
        if ddebug: print (i, distance_pc, magnitude_v, spectrum)
        # It's only meaningful if magnitude_v < 9 because the 
        #   hygxyz.csv is a subset of stars with mag_v < 9.
        if ( distance_pc < crit_distance_pc) and (magnitude_v < crit_magnitude_v) and (spectrum == 'F' or spectrum == 'G' or spectrum == 'K') :
            hygstar_selected.append(hygstar_all[i])
if debug: print ('    the hygstar_selected list contains', len(hygstar_selected), 'stars.\n')
if debug: print ()


# Write the all the selected hyg stars (with 20 pc and < 10 Mag).
write_all_hyg_stars(hygstar_selected, 'hygselected.csv')



######################### The final stage ##########################
# Check if there are at least 8 reference stars in Hipparcos Catalog


# Create the dirctory that will be used to save all the information 
#       of the reference stars for each target
directory = 'reference_stars_of_each_target'
if os.path.exists(directory):
    if debug:
        print ('!! ATTENTION: delete the old "%s" directory after ---\n' % directory )
        print ('         3 seconds')
        time.sleep(1)
        print ('         2 seconds')
        time.sleep(1)
        print ('         1 seconds\n')
        time.sleep(1)
    shutil.rmtree(directory)
    os.makedirs(directory)
else:
    os.makedirs(directory)


if debug: print ('Searching for the reference stars for each body in the hygstar_selected:\n')

# for progress tracking
if debug:
    ten_percent = len(hygstar_selected)/10.0
    n_ten_percent = 1

# all the final target
target_list = []

for i in range(len(hygstar_selected)):
    
    ##### for progress tracking
    if debug:
        if i > int(ten_percent*n_ten_percent):
            print ("    progress: ", n_ten_percent*10, "% have been completed;")
            n_ten_percent = n_ten_percent + 1
    
    ##### Check the reference star for each target. >>> Get the target in the Hipparcos Catalog
    target_once_hyg_id = int(hygstar_selected[i].hyg_id)
    if ddebug: print (i, target_once_hyg_id)

    target_once = hygstar_selected[i]

    # Get the position, parallax, and proper motion of the target 
    target_once_ra24 = target_once.ra_epoch2000
    target_once_ra = float(target_once_ra24) * 15.0
    target_once_dec = float(target_once.dec_epoch2000)
    target_once_pm_ra = float(target_once.pm_ra)
    target_once_pm_dec = float(target_once.pm_dec)
    target_once_magv = target_once.magnitude_v
    target_once_hdname = target_once.hd_name
    target_once_glname = target_once.gliese_name


    fra0dc0 = open('ra0dc0', 'w')
    fra0dc0.write(str(target_once_ra))
    fra0dc0.write('\n')
    fra0dc0.write(str(target_once_dec))
    fra0dc0.close()

    FNULL = open(os.devnull, 'w')
    #subprocess.call("./u4access", "/dev/null")
    subprocess.call(['./u4access'], stdout=FNULL, stderr=subprocess.STDOUT)

    # remove the header of u4.tab
    header_u4 = 9
    with open('u4.tab') as oldfile, open('u4tmp.tab', 'w') as newfile:
        for i, line in enumerate(oldfile):
            if i >= header_u4:
                newfile.write(line)

    os.system('./awk.sh')

    refstar_racDdec_array = np.genfromtxt('u4new.tab', dtype=float)

    n_refstar_racDdec_array, m_refstar_racDdec_array = refstar_racDdec_array.shape

    pm_refstar_array =  np.zeros((n_refstar_racDdec_array,1), dtype=float)

    pm_refstar_array[:,0] = ( refstar_racDdec_array[:,0]**2.0 + refstar_racDdec_array[:,1]**2 )**0.5
    #pm_refstar_array[:,0] = math.sqrt( refstar_racDdec_array[:,0]**2.0 + refstar_racDdec_array[:,1]**2 )

#    pm_refstar_array[:,0] = abs(refstar_array[:,0]) / abs( refstar_array[:,0] - target_once_pm_ra*math.cos(target_once_dec) )
#    pm_refstar_array[:,1] = abs(refstar_array[:,1]) / abs( refstar_array[:,1] - target_once_pm_dec ) 

    pm_target = ( (target_once_pm_ra*math.cos(target_once_dec))**2 + target_once_pm_dec**2 )**0.5

    n_reference = []

    for i in range(n_refstar_racDdec_array):
        if crit_pm:
            if ( (pm_refstar_array[i,0] / pm_target) < crit_pm_ratio):
                n_reference.append(i)
                if ddebug:
                    print (n_reference)
                    print (pm_refstar_array[i,0],pm_target,pm_refstar_array[i,0]/pm_target)
        else:
            n_reference.append(i)
            if ddebug:
                print (n_reference)
                print (pm_refstar_array[i,0],pm_target,pm_refstar_array[i,0]/pm_target)

    if len(n_reference) > 8:
        if ddebug:
            print (n_reference)

        with open('u4tmp.tab') as oldfile, open('u4ref.tab', 'w') as newfile:
            for i, line in enumerate(oldfile):
                for n in range(len(n_reference)):
                    if i == n_reference[n]:
                        newfile.write(line)

        # Save the hip2_id of the target star in the 
        #    './reference_stars_of_each_target/hip2id_(idnumber)' file
        filename_target_once = './reference_stars_of_each_target/hyg_id_%s' % str(target_once_hyg_id)
        file_target_once = open(filename_target_once,'w')
        file_target_once.write("Target star: hyg_id,hipparcos_id,hd_name,hr_name,gliese_name,bayerflamsteed_name,proper_name,ra_epoch2000,dec_epoch2000,distance_pc,pm_ra,pm_dec,rv,magnitude_v,absmagnitude_v,spectrum,color_index,x_epoch2000,y_epoch2000,z_epoch2000,v_x,v_y,v_z\n")
        file_target_once.write('%s\n' % target_once)
        # Save the hip2_id of all the reference stars 
        file_target_once.write("Reference star: from ucac4 Catalog\n")
        file_target_once.write('\n%s\n\n' % n_reference)
        os.system('mv u4ref.tab "%s"_ref ' % filename_target_once )
        #os.system('cat u4.tab >> "%s" ' % filename_target_once )
        #subprocess.call('cat u4.tab >> filename_target_once', shell=True)

        target_list.append(target_once)
        os.system('rm -f u4new.tab')



if debug:
    print ("    progress:  Done.\n")
    print ()

# Save all the targets.
write_all_hyg_stars(target_list, 'targets.csv')
if debug:
    print ( 'Total:  %s targets have been saved in "targets.csv"' % len(target_list) )
    print ( 'The reference stars of each target are saved in the reference_stars_of_each_target directory\n')
    print ( 'May the Targets be with you.\n')




print ('         3 seconds\n')
time.sleep(3)

mult_target = []

targets_num = len(target_list)
for i in range(targets_num):
    one_ra24 = target_list[i].ra_epoch2000
    one_ra = float(one_ra24) * 15.0
    one_dec = float(target_list[i].dec_epoch2000)
    if debug: print (i, one_ra, one_dec)
    for j in range(targets_num):
        if ( i != j):
            another_ra24 = target_list[j].ra_epoch2000
            another_ra = float(another_ra24) * 15.0
            another_dec = float(target_list[j].dec_epoch2000)
            if debug: print (j, another_ra, another_dec)
            if ( ((another_ra - one_ra)**2 + (another_dec - one_dec)**2)**0.5 < 0.44):
                print ( ((another_ra - one_ra)**2 + (another_dec - one_dec)**2)**0.5 )
                #mult_target.append(target_list[j])
                mult_target.append(target_list[i])

# Save all the mult_targets.
write_all_hyg_stars(mult_target, 'mult_target.csv')

