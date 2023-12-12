#!/usr/bin/env bash

NR_NODES=1
TASKS_PER_NODE=24

while getopts ":n:t:" o; do
    case "${o}" in
        n)
            NR_NODES=${OPTARG}
            ;;
        t)
            TASKS_PER_NODE=${OPTARG}
            ;;
    esac
done
shift $((OPTIND-1))

CMD=$(readlink -f $1)
shift

JULIA=$(which julia)

DIR=$(dirname "${BASH_SOURCE[0]}")
JULIA_SETUP=$(readlink -f ${DIR}/slurmSetup.jl)

NR_TASKS=$((NR_NODES * TASKS_PER_NODE))
SCRIPT=__$(basename $CMD)__

cat >$SCRIPT <<EOF
#!/usr/bin/env bash
#
#SBATCH --nodes=$NR_NODES
#SBATCH --ntasks=$NR_TASKS
#SBATCH --ntasks-per-core=1
#
#SBATCH --partition=HIall
#SBATCH --gres=reconExclusive:1
#SBATCH --job-name=$(basename $CMD)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user='olaf.delgado@gmail.com'
#

source ~/.bashrc

$JULIA -L $JULIA_SETUP $CMD $*
EOF

sbatch $SCRIPT

# -- The following is currently useless since our setup script eats stdout
#JOB_ID=$(sbatch --parsable $SCRIPT)
#exec tail -F -n +0 slurm-${JOB_ID}.out
