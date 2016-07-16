#!/bin/sh

# Give data file names as arguments for the entire script

# Set slurm log directory (make sure ~ is expanded)
slurm_log_dir=`echo ~/simulation_code_bin/_slurm_files/`

# Set number of cores to use
n_cores=4

# Specify the real memory required per node in MegaBytes.
mem_per_node=24000

# Maximum runtime in minutes
run_time=7200

# SLURM partitions
# https://rc.fas.harvard.edu/resources/running-jobs/#SLURM_partitions
partition="serial_requeue"


# Define a function with two arguments
RunScript() {
    # Argument 1 is data file name (invoke from the script directory)
    data_file=$1
    # Argument 2 is R script name
    r_script=$2


    # Generate slurm log file names (remove path)
    slurm_out_txt=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}.out.txt
    slurm_err_txt=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}.err.txt

    # Generate script name
    script_name=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}.slurm


    # This part creates the script
    # This flushes the file if it exists
    echo "#!/bin/bash" > ${script_name}

    echo "#SBATCH -n "${n_cores}" # Number of cores requested" >> ${script_name}
    echo "#SBATCH -N 1 # Ensure that all cores are on one machine" >> ${script_name}
    echo "#SBATCH -t "${run_time}" # Runtime in minutes" >> ${script_name}
    echo "#SBATCH -p "${partition}" # Partition to submit to" >> ${script_name}
    echo "#SBATCH --mem="${mem_per_node}" # real memory required per node" >> ${script_name}
    echo "#SBATCH -o "${slurm_log_dir}${slurm_out_txt}" # Standard out goes to this file" >> ${script_name}
    echo "#SBATCH -e "${slurm_log_dir}${slurm_err_txt}" # Standard err goes to this file" >> ${script_name}

    echo "# Load R" >> ${script_name}
    echo "# module load R_packages/3.2.2-fasrc01" >> ${script_name}
    echo "module load centos6/R-3.1.2" >> ${script_name}

    echo "# Configure R local package directory AFTER R has been loaded" >> ${script_name}
    echo "export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER" >> ${script_name}

    echo "# Invoke simulation runner with file name" >> ${script_name}
    echo "Rscript" ${r_script} ${data_file} >> ${script_name}

    # Show script
    echo ""
    echo "Running this script"
    cat ${script_name}
    echo ""

    # This part runs the script (data file name is coded)
    echo sbatch ${script_name}
    sbatch ${script_name}

    # This part removes the script
    echo rm ${script_name}
    rm ${script_name}
}


for file in $@
do
    RunScript ${file} 02.RunSimulation.R
done
