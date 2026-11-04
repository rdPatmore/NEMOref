#!/bin/bash
#SBATCH --job-name=nemo00
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --account=n02-ncas
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos
# Created by: mkslurm -S 1 -s 1 -m 1 -C  1 -c  1 -t 00:10:00 -a n01 -j nemo_test
module -s restore /work/n01/shared/acc/n01_modules/ucx_env
export OMP_NUM_THREADS=1
#
cat > myscript_wrapper2.sh << EOFB
#!/bin/ksh
#
set -A map ./xios_server.exe ./nemo
exec_map=( 0 1 )
#
exec \${map[\${exec_map[\$SLURM_PROCID]}]} 
##
EOFB
chmod u+x ./myscript_wrapper2.sh
#
srun --mem-bind=local --cpu-bind=v,map_cpu:00,0x1, ./myscript_wrapper2.sh
