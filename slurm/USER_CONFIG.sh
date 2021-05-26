
if [ $USER == littlem ]; then
        SCRIPTS_FOLDER="/cs/labs/ravehb/littlem/npctransport_kinetic/"

elif [ $USER == ravehb ]; then
    SCRIPTS_FOLDER="/cs/labs/ravehb/ravehb/Projects/npctransport_kinetic/"
else
    echo "WARNING: Unknown user - using ravehb config"
    SCRIPTS_FOLDER="/cs/labs/ravehb/ravehb/Projects/npctransport_kinetic/"
fi

SLURM_FOLDER="$SCRIPTS_FOLDER/slurm"
