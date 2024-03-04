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
print("Hello & good luck optimizing the cloud profiles! The haze parameters are already set to 4 ppm for a titan like haze with a particle size of 0.6 microns (these values are determined in optimize plots). As a starting point, consider that NH4SH clouds exist between approximately 120-200 Kelvin.")

folder='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/'
# direction to the default jupiter configuration file from PSG
# the jupiter compatible file also includes the type of output data that we will want (i/f albedo) and spectral range 0.2-3 micron (among other variations)
jupiter_base_cfg= open(folder+"psg_cfg_2ppm_0.05haze.txt", "r");

# reads/ stores each line of file in basecfg (as a list)
basecfg= jupiter_base_cfg.readlines();
#close the file everytime
jupiter_base_cfg.close();

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

print("MOLECULES in this template file:")
for value in molecules:
    print(value)

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

#                       TO RUN THE FUNCTION ON PSG AND GET A RETURN PLOT

# -------------------------------------------------------------------------------------------------------------
diff_vals=[]
def change_cloud_height_psgspectrum(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # consider this a scale factor for the cloud abundance, it will be automatically set to 1 unless specified in function
    # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
    print("Running cloud function...")
    for value in range(len(molecules)):
        molecule_index=0
        if molecules[value][0:2]== 'C2':
            if molecules[value][0:4]== absorber[0:4]:
                molecule_index= value
                break
        if molecules[value][0:3]== absorber[0:3]:
            molecule_index= value
            #stop the loop once this requirement is satisfied
            break
    
    # This is where the true data from the ~original~ config file is stored 
    store_data_molecules_original = np.zeros(shape=(num_layers_int, len(molecules)+2))
    loc_layers=0
    for line in basecfg:
        if "<ATMOSPHERE-LAYER-" in line:
            abun_1= line.split(',')
            # adjust the formatting from the original config file
            abun_1[0]= abun_1[0].split('>')[1]
            #store_data_molecules[:,]
            loc_layers= loc_layers+1
            # this line stores all of the data into an array that can be adjusted to change the height of the cloud deck
            store_data_molecules_original[loc_layers-1]= abun_1
    # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
            #the clouds of a given absorber 
            # need to add 2 to the indexing bc P&T are the first 2 values
    max_value_multiply= np.max(store_data_molecules_original[:,molecule_index+2])

    # copy the data frame so this version can be modified and compared to the original one
    store_data_molecules_newvals= np.zeros(shape=(num_layers_int, len(molecules)+2))
    for value in np.arange(0,len(store_data_molecules_original),1):
        store_data_molecules_newvals[value]= store_data_molecules_original[value]
   

    # Make the config file:
    file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_absorber'+str(absorber)+'_min_'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'.txt'
    f = open(file_name, 'w')
    # copy the original jupiter file to vary the values before saving
    new_file= basecfg.copy() 
    
    for line_val in np.arange(min_level, max_level):
        # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
        # the clouds of a given absorber
        new_val= (store_data_molecules_original[line_val,(molecule_index+2)]+ scale_abun*max_value_multiply)
        store_data_molecules_newvals[line_val,(molecule_index+2)]= new_val
    
    for value in np.arange(min_level, max_level,1): 
        new_file[value+loc_layer1-1]= "<ATMOSPHERE-LAYER-" + str(value) + ">"
        line_values = ','.join(map(str, store_data_molecules_newvals[value]))
        new_file[value+loc_layer1-1] += line_values + '\n'

    for line in new_file:
        f.write(str(line))
    f.close() 

    print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2])
    print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2])

    # RUN THROUGH PSG
    print("Sending amazing science to PSG...")
    print("  ")

    run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
    os.system(run_file)
    print("Good news, it worked!")

    #  ~~~~~~~~~~~~~~ Plotting + Optimize
    # Read in data files
    #jupiter_published= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/madden_spectra/CatalogofSolarSystemObjects/Albedos/Jupiter_Lundock080507_Albedo.txt')
    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
    output_file_NEW_psg=np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt')
    # Read in Karkoschka data
    day_night_jup= np.genfromtxt('../jupiter_spectra/coulterbarnesfortney_fig5_jup_dbf5.txt')
    karkoschka_fromgv= np.genfromtxt('../jupiter_spectra/KarkoschkaAlbedoJupiter.txt')
    # The data for each of the spectral data types (day/night, cloudy/clear) come from day_night_jup. need to sort:
    dayside_cloudy= np.array(day_night_jup[0:330])

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
    plt.figure(dpi=200)
    plt.plot((karkoschka_fromgv[:,0])/1000, karkoschka_fromgv[:,3], label='Karkoschka', linewidth=1, color='royalblue')
    plt.plot(dayside_cloudy[:,1],dayside_cloudy[:,3], label='Dayside Cloudy CBF (MAIN REFERENCE)', linewidth=0.7, color='dodgerblue')
    plt.plot(output_file_NEW_psg[:,0], output_file_NEW_psg[:,1], label='PSG NEW Model; RMSE:{}'.format(new_psg_vs_cbf), color= 'limegreen')
    plt.plot(psg_default_data_jupiter[:,0], psg_default_data_jupiter[:,1], label='PSG Default; RMSE:{}'.format(old_psg_vs_cbf), linewidth=0.7, linestyle='--', color='violet')
    plt.grid()
    plt.axvline(x=.678, label='Peak I/F occurs at 0.678 microns', color='sandybrown', linestyle='-.', linewidth=0.5)
    plt.ylim(0,0.6)
    plt.xlim(0,4.8)
    plt.title('Zoomed on Karkoschka Data; {} Fit; {} Cloud at layers {}:{}, Abundance: x{}'.format(title, absorber, min_level, max_level, scale_abun))
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
    
    #plt.savefig('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_plots/'+'Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_PLOT.png');
    
    return print("The change has made the fit {}".format(title))

#change_cloud_height_psgspectrum(absorber='NH4SH', scale_abun=4, min_level=5, max_level=10)