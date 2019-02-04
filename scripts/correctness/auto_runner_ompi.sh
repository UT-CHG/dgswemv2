#!/bin/bash

export OMP_NUM_THREADS=4

declare -a meshsizes=(1)
declare -a polyorders=(1)

for i in ${meshsizes[@]}; do 
	for j in ${polyorders[@]}; do
	        #make mesh
		cd ~/../../workspace/dgswemv2/build/mesh_generators/
		sed -i "7s/.*/num_x_subdivisions: $i/" rectangle_manufactured
		sed -i "8s/.*/num_y_subdivisions: $i/" rectangle_manufactured
		./rectangular_mesh_generator rectangle_manufactured 
		mv rectangular_mesh.14 ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution
	
		#simulate using mesh
		cd ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution/
        	sed -i "19s/.*/polynomial_order: $j/" dgswemv2_input.15
		
		#parallelize mesh
		 cd ~/../../workspace/dgswemv2/build/partitioner/ 
		./partitioner ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution/dgswemv2_input.15 4 1
		
		cd ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution
		echo "Meshsize = $i; Polynomial Order = $j " >> ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt
        	./RKDG_MANUFACTURED_SOLUTION_OMPI dgswemv2_input_parallelized.15 |tee -a ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt
		echo ' ' >> ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt

		rm rectangular_mesh*
	done 
done
 
