./Allclean
blockMesh
decomposePar
sbatch jobAzure.sh
#NPROCS=$(awk '/numberOfSubdomains/{gsub(";","",$2); print $2}' system/decomposeParDict)
#/shared/home/azureuser/lib/openfoam/mto/ThirdParty-6/platforms/linux64Gcc/openmpi-1.10.7/bin/mpirun -np "$NPROCS" MTO_ThermalFluid -parallel > log &
