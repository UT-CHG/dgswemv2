#!/bin/bash

declare -a meshsizes=(4 8 16 32 64)
declare -a polyorders=(1 2 3 4 5)
declare -a timesteps=(4 2 1 0.5 0.25) 

for i in ${meshsizes[@]}; do 
    for j in ${polyorders[@]}; do
	for k in ${timesteps[@]}; do 
	        #make mesh
	       	cd ~/../../workspace/dgswemv2/build/mesh_generators/
	        sed -i "7s/.*/num_x_subdivisions: $i/" rectangle_manufactured
	        sed -i "8s/.*/num_y_subdivisions: $i/" rectangle_manufactured
	        ./rectangular_mesh_generator rectangle_manufactured 
        	mv rectangular_mesh.14 ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution
        	#simulate using mesh
        	cd ~/../../workspace/dgswemv2/build/examples/rkdg_swe_manufactured_solution/
		sed -i "15s/.*/dt: $k/" dgswemv2_input.15 
        	sed -i "19s/.*/polynomial_order: $j/" dgswemv2_input.15
	#	echo "Meshsize = $i; Polynomial Order = $j; Timestep = $k: " >> ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt
        	./RKDG_MANUFACTURED_SOLUTION_SERIAL dgswemv2_input.15# |tee -a ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt
	#	echo "\n" >> ~/../../workspace/dgswemv2/build/examples/output/timenerror.txt
	done 
    done 
done 
