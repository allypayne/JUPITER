import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------------------------------------------

#                              READING IN DATA FOR ALL FUNCTIONS/ SYTNAX

# -------------------------------------------------------------------------------------------------------------

folder='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/'
# direction to the default jupiter configuration file from PSG
# the jupiter compatible file also includes the type of output data that we will want (i/f albedo) and spectral range 0.2-3 micron (among other variations)
jupiter_base_cfg= open(folder+"psg_cfg_2ppm_0.05haze.txt", "r");

# reads/ stores each line of file in basecfg (as a list)
basecfg= jupiter_base_cfg.readlines();
#close the file everytime
jupiter_base_cfg.close();

# Read in Karkoschka data
day_night_jup= np.genfromtxt('../jupiter_spectra/coulterbarnesfortney_fig5_jup_dbf5.txt')
karkoschka_fromgv= np.genfromtxt('../jupiter_spectra/KarkoschkaAlbedoJupiter.txt')
# The data for each of the spectral data types (day/night, cloudy/clear) come from day_night_jup. need to sort:
dayside_cloudy= np.array(day_night_jup[0:330])



## (Here we find "loc_molecules"- line where the base_config lists the molecules that) are considered in jupiters atmosphere (use this index to print the molecules)
# BASE CONFIG: "<ATMOSPHERE-LAYERS-MOLECULES>H2O,CO2,O3..."
loc_molecules=0
for line in basecfg:
    if "ATMOSPHERE-LAYERS-MOLECULE" in line:  
        break
    loc_molecules= loc_molecules+1

molecules=[]
for line in basecfg:
    if '<ATMOSPHERE-LAYERS-MOLECUL' in line:
        molecule= line.split(',')
        molecule[0]= molecule[0][-4:]
        for value in molecule:
            molecules.append(value)

print("Welcome to the Hazes Function!")

# keep track of the # of layers here
for line in basecfg:
    if "<ATMOSPHERE-LAYERS>" in line:
        num_layers= line.split('>')[1]
        num_layers_int= int(num_layers)

loc_layer1=0
for line in basecfg:
    if "<ATMOSPHERE-LAYER-1" in line:  
        break
    loc_layer1= loc_layer1+1
    
coulter_barnes_fortney= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/coulterbarnesfortney_fig4_jup.txt')
    
phase_zero_wavelength=[]
phase_zero_albedo=[]
for row in coulter_barnes_fortney:
    if row[0]<2:
        phase_zero_wavelength.append(row[1])
        phase_zero_albedo.append(row[3])
    
# Molecules included: 'C2H2', 'C2H4', 'C2H6', 'CH4', 'H2O', 'NH3', 'WaterIce', 'NH4SH', 'Ammonia'

# -------------------------------------------------------------------------------------------------------------

#            TO RUN A RANGE OF FUNCTIONS ON PSG AND FIND THE PARAMETERS THAT IMPROVE THE FIT THE MOST

# -------------------------------------------------------------------------------------------------------------

diff_vals=[]
def vary_particle_params(abund_min, abund_max, min_size, max_size):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # this is the range of abundance scale factors we want to consider (units of ppm)
    for value_1 in np.arange(abund_min,abund_max,1):
        # this is the range of particle sizes we want to consider (in units of micron)
        for value_2 in np.arange(min_size,max_size,0.01):

            # Make the config file:
            file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_ppmabundance_'+str(value_1)+"_micronsize_"+str(value_2)+'.txt'
            f = open(file_name, 'w')
            # copy the original jupiter file to vary the values before saving
            new_file= basecfg.copy() 

            # Define the locations of the lines that I want to change
            loc_aabun=0
            for line in basecfg:
                if "<ATMOSPHERE-AABUN" in line:  
                    break
                loc_aabun= loc_aabun+1
            print(loc_aabun)

            loc_asize=0
            for line in basecfg:
                if "<ATMOSPHERE-ASIZE" in line:  
                    break
                loc_asize= loc_asize+1
            print(loc_asize)

            for line in basecfg:
                if "<ATMOSPHERE-AABUN" in line:
                    line_new_1= '<ATMOSPHERE-AABUN>1,1,'+str(value_1)+'\n'
                    #print(line_new_1)

                if "<ATMOSPHERE-ASIZE" in line:
                    line_new_2= "<ATMOSPHERE-ASIZE>0.1,0.01,"+str(value_2)+'\n'
                    #print(line_new_2)


            new_file[loc_aabun]= line_new_1
            new_file[loc_asize]= line_new_2

            for line in new_file:
                f.write(str(line))
            f.close() 

            # RUN THROUGH PSG
            print("Sending amazing science to PSG...")
            print("  ")

            run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./output_psg_spectra/'+'Jupiter_ppmabundance_'+str(value_1)+"micronsize"+str(value_2)+'_OUT.txt'
            os.system(run_file)
            print("It worked!")
            
        #  ~~~~~~~~~~~~~~ Plotting + Optimize
            # Read in data files
            #jupiter_published= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/madden_spectra/CatalogofSolarSystemObjects/Albedos/Jupiter_Lundock080507_Albedo.txt')
            psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
            output_file_NEW_psg=np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_ppmabundance_'+str(value_1)+"micronsize"+str(value_2)+'_OUT.txt')

            # the spectrum that we want to match is CBF Cloudy
            jupiter_ideal_x= dayside_cloudy[:,1]
            jupiter_ideal_y= dayside_cloudy[:,3]

            psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
            x_vals_default= psg_default_data_jupiter[:,0]
            y_vals_default= psg_default_data_jupiter[:,1]

            # Cubic Spline functions
            f_cbf_initial_spectrum = CubicSpline(x_vals_default,y_vals_default, bc_type='natural')
            f_cbf_new_data_spectrum = CubicSpline(output_file_NEW_psg[:,0],output_file_NEW_psg[:,1], bc_type='natural')

            new_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_new_data_spectrum(jupiter_ideal_x))**2)/(len(jupiter_ideal_y)-2))
            old_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_initial_spectrum(jupiter_ideal_x))**2)/(len(jupiter_ideal_y)-2))
            delta_rmse= old_psg_vs_cbf-new_psg_vs_cbf

            diff_vals.append(delta_rmse)

            if new_psg_vs_cbf<old_psg_vs_cbf:
                print("This change improved the fit!")
                print("  ")
                title='BETTER'
                # COULD ADD LINE HERE TO MAKE THE BASE CFG FOR THE NEXT RUN=OUTPUT OF THIS ONE (IF IMPROVED)
            else:
                print("The new fit is worse than the default data.")
                title='WORSE'
            print("New data vs ideal:",new_psg_vs_cbf)
            print("Default data vs ideal:", old_psg_vs_cbf)

            ###
            # ~~~~~~~~~~~~~~~~~ The small/ zoomed plot is here ~~~~~~~~~~~~~~~~
            ###

            plt.plot((karkoschka_fromgv[:,0])/1000, karkoschka_fromgv[:,3], label='Karkoschka', linewidth=1, color='royalblue')
            plt.plot(dayside_cloudy[:,1],dayside_cloudy[:,3], label='Dayside Cloudy CBF (MAIN REFERENCE)', linewidth=0.7, color='dodgerblue')
            plt.plot(output_file_NEW_psg[:,0], output_file_NEW_psg[:,1], label='PSG NEW Model; RMSE:{}'.format(new_psg_vs_cbf), color= 'limegreen')
            plt.plot(psg_default_data_jupiter[:,0], psg_default_data_jupiter[:,1], label='PSG Default; RMSE:{}'.format(old_psg_vs_cbf), linewidth=0.7, linestyle='--', color='violet')
            plt.grid()
            plt.axvline(x=.678, label='Peak I/F occurs at 0.678 microns', color='sandybrown', linestyle='-.', linewidth=0.5)
            plt.ylim(0,0.6)
            plt.xlim(0,4.8)
            plt.title('Zoomed on Karkoschka Data; {} Fit; PSG Haze Abundance {}, Particle Size: {}'.format(title, value_1, value_2))
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel('I/F')
            plt.xticks(np.arange(0,4.8, 0.2))
            plt.grid(True, linewidth=0.4, linestyle='--') 
            plt.legend(fontsize=7)
            plt.xlim(0.3,1)
            plt.show()
            
            ###
            # The larger plot is here 
            ###
            plt.figure(figsize=(9,5), dpi=200)
            plt.plot((karkoschka_fromgv[:,0])/1000, karkoschka_fromgv[:,3], label='Karkoschka', linewidth=1, color='royalblue')
            plt.plot(dayside_cloudy[:,1],dayside_cloudy[:,3], label='Dayside Cloudy CBF (MAIN REFERENCE)', linewidth=0.7, color='dodgerblue')
            plt.plot(output_file_NEW_psg[:,0], output_file_NEW_psg[:,1], label='PSG NEW Model; RMSE:{}'.format(new_psg_vs_cbf), color= 'limegreen')
            plt.plot(psg_default_data_jupiter[:,0], psg_default_data_jupiter[:,1], label='PSG Default; RMSE:{}'.format(old_psg_vs_cbf), linewidth=0.7, linestyle='--', color='violet')
            plt.ylim(0,0.7)
            plt.xlim(0,4.8)
            plt.title('Jupiter Spectra')
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel('I/F')
            plt.xticks(np.arange(0,4.8, 0.5))
            plt.grid(True, linewidth=0.4, linestyle='--') 
            plt.legend(fontsize=7)
            plt.text(2, 0.33, 'Features here are more impacted by larger particle sizes', fontsize=7, color='blue')
            plt.text(0.1, 0.57, 'This side is more impacted by', fontsize=7, color='blue')
            plt.text(0.1, 0.55, 'smaller particles', fontsize=7, color='blue')
            plt.show();

    return diff_vals


#vary_particle_params(abund_min=1, abund_max=4, min_size=0.01, max_size=.11)