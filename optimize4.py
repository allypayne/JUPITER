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
jupiter_base_cfg= open(folder+"psg_cfg_5micron_compat_test.txt", "r");

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

#                       IF YOU JUST WANT THE INPUT/ CONFIG FILES RUN THIS FUNCTION

# -------------------------------------------------------------------------------------------------------------
def change_cloud_height_getfile(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # consider this a scale factor for the cloud abundance, it will be automatically set to 1 unless specified in function
    # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
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
    #np.savetxt('original_molecules.csv', rounded_data, delimiter=',')

    # need to add 2 to the indexing bc P&T are the first 2 values
    max_value_multiply= np.max(store_data_molecules_original[:,molecule_index+2])
    # copy the data frame so this version can be modified and compared to the original one
    store_data_molecules_newvals= np.zeros(shape=(num_layers_int, len(molecules)+2))
    for value in np.arange(0,len(store_data_molecules_original),1):
        store_data_molecules_newvals[value]= store_data_molecules_original[value]
   
    # Make the config file:
    file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'.txt'
    f = open(file_name, 'w')
    # copy the original jupiter file to vary the values before saving
    new_file= basecfg.copy() 
    
    for line_val in np.arange(min_level, max_level):
        # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
        # the clouds of a given absorber
        new_val= (store_data_molecules_original[line_val,(molecule_index+2)]+ scale_abun*max_value_multiply)
        store_data_molecules_newvals[line_val,(molecule_index+2)]= new_val

    # New vals has length 30 (This corresponds to the 30 changed lines; does not have data for every layer)
    # for line in store_data_molecules_newvals:
    #     print("NEWVALS LINES:", line)
    
    for value in np.arange(min_level, max_level,1): 
        new_file[value+loc_layer1-1]= "<ATMOSPHERE-LAYER-" + str(value) + ">"
        line_values = ','.join(map(str, store_data_molecules_newvals[value]))
        new_file[value+loc_layer1-1] += line_values + '\n'

    for line in new_file:
        f.write(str(line))
    f.close() 

    return print(r"The atmospheric abundance values for {} have been changed...".format(molecules[molecule_index])), print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2]), print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2]), print("The file is saved as", file_name)

# Run line below to try this function
#change_cloud_height_getfile(absorber='Ammonia', scale_abun=20, min_level=30, max_level=45)
# -------------------------------------------------------------------------------------------------------------

#                   TO RUN THE FUNCTION ON PSG AND SEE IF YOUR CHANGES IMPROVED OR WORSEN THE FIT

# -------------------------------------------------------------------------------------------------------------

def change_cloud_height_psgspectrum(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
    print("Running cloud function...")
    for value in range(len(molecules)):
        molecule_index=0
        if molecules[value][0:1]== 'C2':
            if molecules[value][0:4]== absorber[0:4]:
                molecule_index= value
                break    
        if molecules[value][0:2]== absorber[0:2]:
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
    file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'.txt'
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
    #curl --data-urlencode file@config.txt https://psg.gsfc.nasa.gov/api.php
    run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
    #os.system(run_file)
    print("Good news, it worked!")
    
  #  ~~~~~~~~~~~~~~~~~ Plotting + Optimize
    # Read in data files
    jupiter_published= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/madden_spectra/CatalogofSolarSystemObjects/Albedos/Jupiter_Lundock080507_Albedo.txt')
    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
    output_file_psg='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
    psg_NEW_data= np.genfromtxt(output_file_psg)
    new_fit_psg_x= psg_NEW_data[:,0]
    new_fit_psg_y= psg_NEW_data[:,1]

    jupiter_ideal_x= phase_zero_wavelength
    jupiter_ideal_y= phase_zero_albedo

    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
    x_vals_default= psg_default_data_jupiter[:,0]
    y_vals_default= psg_default_data_jupiter[:,1]

    # Cubic Spline functions
    f_cbf_initial_spectrum = CubicSpline(x_vals_default,y_vals_default, bc_type='natural')
    f_cbf_new_data_spectrum = CubicSpline(new_fit_psg_x,new_fit_psg_y, bc_type='natural')

    new_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_new_data_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
    old_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_initial_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
    delta_rmse= old_psg_vs_cbf-new_psg_vs_cbf
    if new_psg_vs_cbf<old_psg_vs_cbf:
        print("This change improved the fit! (:")
        print("  ")
        title='BETTER'
        # COULD ADD LINE HERE TO MAKE THE BASE CFG FOR THE NEXT RUN=OUTPUT OF THIS ONE (IF IMPROVED)
    else:
        print("The new fit is worse than the default data.")
        title='WORSE'
    print("New data vs ideal:",new_psg_vs_cbf)
    print("Default data vs ideal:", old_psg_vs_cbf)

    plt.figure(dpi=250)
    plt.plot(psg_default_data_jupiter[:,0],psg_default_data_jupiter[:,1], label='PSG DEFAULT', linestyle='--', color='mediumturquoise')
    plt.plot(jupiter_ideal_x,jupiter_ideal_y, label='CBF (Ideal) Spectrum')
    plt.plot(new_fit_psg_x,new_fit_psg_y, label='New fit', linewidth=0.4, color='forestgreen')    
    plt.grid()
    plt.legend()
    plt.ylabel('I/F')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.suptitle('{} ; $\Delta$ RMSE: {}'.format(title, delta_rmse))
    plt.title('Adjusted Jupiter Params: {} clouds at layers {}-{}, scaled by {}x'.format(molecules[molecule_index],str(min_level),str(max_level),str(scale_abun)))
    plt.ylim(0,1)
    plt.show()

#change_cloud_height_psgspectrum('CH4', scale_abun=5, min_level=30, max_level=45)

# -------------------------------------------------------------------------------------------------------------

#            TO RUN A RANGE OF FUNCTIONS ON PSG AND FIND THE PARAMETERS THAT IMPROVE THE FIT THE MOST

# -------------------------------------------------------------------------------------------------------------


def optimize_cloud_height_psgspectrum(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    diff_vals=[]
    # TREAT this as min level (value_1)
    for value_1 in range(min_level,max_level):
        # TREAT this as max level (value_2)
        for value_2 in range(min_level,max_level):
            # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
            print("Running cloud function...")
            for value in range(len(molecules)):
                molecule_index=0
                if molecules[value][0:1]== 'C2':
                    if molecules[value][0:4]== absorber[0:4]:
                        molecule_index= value
                        break    
                if molecules[value][0:2]== absorber[0:2]:
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
            file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_min'+str(value_1)+'_max'+str(value_2)+'_scale'+str(scale_abun)+'.txt'
            f = open(file_name, 'w')
            # copy the original jupiter file to vary the values before saving
            new_file= basecfg.copy() 
            
            for line_val in np.arange(value_1, value_2):
                # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
                # the clouds of a given absorber
                new_val= (store_data_molecules_original[line_val,(molecule_index+2)]+ scale_abun*max_value_multiply)
                store_data_molecules_newvals[line_val,(molecule_index+2)]= new_val
            
            for value in np.arange(value_1, value_2,1): 
                new_file[value+loc_layer1-1]= "<ATMOSPHERE-LAYER-" + str(value) + ">"
                line_values = ','.join(map(str, store_data_molecules_newvals[value]))
                new_file[value+loc_layer1-1] += line_values + '\n'

            for line in new_file:
                f.write(str(line))
            f.close() 

            #print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2])
            #print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2])

            # RUN THROUGH PSG
            print("Sending amazing science to PSG...")
            print("  ")
            run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(value_1)+'_max'+str(value_2)+'_scale'+str(scale_abun)+'_OUT.txt'
            os.system(run_file)
            print("It worked!")
            
        #  ~~~~~~~~~~~~~~ Plotting + Optimize
            # Read in data files
            jupiter_published= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/madden_spectra/CatalogofSolarSystemObjects/Albedos/Jupiter_Lundock080507_Albedo.txt')
            psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
            output_file_psg='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(value_1)+'_max'+str(value_2)+'_scale'+str(scale_abun)+'_OUT.txt'
            psg_NEW_data= np.genfromtxt(output_file_psg)
            new_fit_psg_x= psg_NEW_data[:,0]
            new_fit_psg_y= psg_NEW_data[:,1]

            jupiter_ideal_x= phase_zero_wavelength
            jupiter_ideal_y= phase_zero_albedo

            psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
            x_vals_default= psg_default_data_jupiter[:,0]
            y_vals_default= psg_default_data_jupiter[:,1]

            # Cubic Spline functions
            f_cbf_initial_spectrum = CubicSpline(x_vals_default,y_vals_default, bc_type='natural')
            f_cbf_new_data_spectrum = CubicSpline(new_fit_psg_x,new_fit_psg_y, bc_type='natural')

            new_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_new_data_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
            old_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_initial_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
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

            plt.figure(dpi=250)
            plt.plot(psg_default_data_jupiter[:,0],psg_default_data_jupiter[:,1], label='PSG DEFAULT', linestyle='--', color='mediumturquoise')
            plt.plot(jupiter_ideal_x,jupiter_ideal_y, label='CBF (Ideal) Spectrum')
            plt.plot(new_fit_psg_x,new_fit_psg_y, label='New fit', linewidth=0.4, color='forestgreen')    
            plt.grid()
            plt.legend()
            plt.ylabel('I/F')
            plt.xlabel('Wavelength ($\mu$m)')
            plt.suptitle('{} ; $\Delta$ RMSE: {}'.format(title, delta_rmse))
            plt.title('Adjusted Jupiter Params: {} clouds at layers {}-{}, scaled by {}x'.format(molecules[molecule_index],str(value_1),str(value_2),str(scale_abun)))
            plt.ylim(0,1)
            plt.show()

    return diff_vals

# To test this function:
optimize_cloud_height_psgspectrum('Ammonia', scale_abun=5, min_level=30, max_level=40)


# Goal of optimize4: improve the for loop to more effectively find solutions.
# have it take the data from the old run that was good and save it to the new cfg file?
# create a for loop that will vary the particle size particularly