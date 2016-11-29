    # Nom arbitraire du travail LoadLeveler
    # @ job_name = run_test_MGROMS
    # Fichier de sortie standard du travail
    # @ output   = $(job_name).$(jobid)
    # Fichier de sortie d'erreur du travail
    # @ error= $(job_name).$(jobid)
    # Type de travail
    # @ job_type = mpich
    # Nombre de processus MPI demandes
    # @ total_tasks = 4
    # Permet le passage de total_tasks a mpirun via NB_TASKS
    # @ environment = NB_TASKS=$(total_tasks)
    # Temps du job hh:mm:ss (10mn ici)
    # @ wall_clock_limit = 00:10:00
    # @ queue

    export TESTDIR='test_testseamount_2x2_v0'
    export SRCDIR='/linkhome/rech/dgw/rdgw004/MGROMS/mgroms-0.3.8/src'
     
    # Recommandation : compilez et exécutez vos codes sous un même environnement Intel.
    # Donc, si necessaire, utilisez la commande module pour charger l'environnement approprie.
    # Par exemple, si votre code est compile avec Intel/2016.2, de-commenter la ligne suivante :
    # module load intel/2016.2

    module load intel/2017.0 netcdf/seq/4.3.3.1 

    export LD_LIBRARY_PATH=/smplocal/pub/NetCDF/4.3.3.1/seq/lib:$LD_LIBRARY_PATH
     
    # Pour avoir l'echo des commandes
    set -x

    # Repertoire temporaire de travail
    cd $WORKDIR

    mkdir ${TESTDIR}
    cd ${TESTDIR}

    # La variable LOADL_STEP_INITDIR est automatiquement positionnee par
    # LoadLeveler au repertoire dans lequel on tape la commande llsubmit
    cp $SRCDIR/testseamount ./.
    cp $LOADL_STEP_INITDIR/nh_namelist ./.

    ls -rtl
    date

    # Execution d'un programme parallèle MPI.
    time mpirun -np $NB_TASKS ./testseamount

    date
    ls -rtl

