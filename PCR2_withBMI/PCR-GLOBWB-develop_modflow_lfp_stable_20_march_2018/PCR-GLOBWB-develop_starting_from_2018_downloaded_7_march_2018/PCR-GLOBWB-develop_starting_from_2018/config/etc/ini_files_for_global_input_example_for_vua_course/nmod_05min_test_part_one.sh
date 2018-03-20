#!/bin/bash
#SBATCH -N 1
#SBATCH -t 119:59:00 
#SBATCH -p normal
#SBATCH --constraint=haswell

# provide the job name
#SBATCH -J test-p1-edwinvua


# mail alert at start, end and abortion of execution                                                                                            
#SBATCH --mail-type=ALL                                                                                                                         
# send mail to this address                                                                                                                     
#SBATCH --mail-user=edwinkost@gmail.com   


cd /home/edwinvua/github/UU-Hydro/PCR-GLOBWB_model/model/
python parallel_pcrglobwb_runner.py /home/edwinvua/github/edwinkost/PCR-GLOBWB/config/ini_files_for_global_input_example_for_vua_course/setup_05min_non-natural_test_with_edwinvua_part_one.ini
