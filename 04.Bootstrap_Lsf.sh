#!/bin/sh

# Give data file names as arguments for the entire script
# Run in the script directory

# Set lsf log directory (make sure ~ is expanded)
lsf_log_dir=`echo ~/simulation_code_bin/_lsf_files/`

# Set number of cores to use
n_cores=8

# Set maximum time
max_time="12:00"

# Set Orchestra queue to use
queue="mcore"


# Define a function with two arguments
RunScript() {
    # Argument 1 is R script name
    r_script=$1
    # Argument 2 is data file name (invoke from the script directory)
    data_file=$2
    # Argument 3 is which part of data to work on
    parti=$3

    # generate a part string
    parti_str=`printf "%0*d\n" 2 ${parti}`

    # Generate lsf log file names (remove path)
    lsf_out_txt=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}_${parti_str}.out.txt
    lsf_err_txt=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}_${parti_str}.err.txt

    # Generate script name
    script_name=`echo ${data_file%\.*} | sed -e "s/.*\///g"`.${r_script%\.*}_${parti_str}.lsf


    # This part creates the script
    # This flushes the file if it exists
    echo "#!/bin/bash" > ${script_name}

    echo "#BSUB -n "${n_cores}" # Number of cores requested" >> ${script_name}
    echo "#BSUB -W "${max_time}" # Job maximum execution time" >> ${script_name}
    echo "#BSUB -J "${script_name}" # Job name" >> ${script_name}
    echo "#BSUB -o "${lsf_log_dir}${lsf_out_txt}" # Standard out goes to this file" >> ${script_name}
    echo "#BSUB -e "${lsf_log_dir}${lsf_err_txt}" # Standard err goes to this file" >> ${script_name}

    # https://wiki.med.harvard.edu/Orchestra/ChoosingAQueue
    # The mcore queue has higher priority than the short or long queues,
    # and jobs in it can't be suspended
    echo "#BSUB -q "${queue}" # Submit to this queue" >> ${script_name}

    # Do not e-mail. Too many messages.
    # echo "#BSUB -N # Also e-mail to the address specified in ~/.forward file" >> ${script_name}
    echo "export OMP_NUM_THREADS="${n_cores} >> ${script_name}

    echo "# Load R" >> ${script_name}
    echo "module load stats/R/3.2.1" >> ${script_name}

    echo "# Configure R local package directory AFTER R has been loaded" >> ${script_name}
    echo "export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER" >> ${script_name}

    echo "# Invoke simulation runner with file name" >> ${script_name}
    echo "Rscript" ${r_script} ${data_file} ${parti} >> ${script_name}

    # Show script
    echo ""
    echo "Running this script"
    cat ${script_name}
    echo ""

    # This part runs the script (data file name is coded)
    bsub < ${script_name}

    # This part removes the script
    echo rm ${script_name}
    rm ${script_name}
}


for file in $@
do
    for i in {1..10}
    do
        # Run i-th part
        RunScript Bootstrap.R ${file} ${i}
    done
done
