#!/bin/bash


cd tp/TP_Poisson_C_students_2022



cp ambre.mk ./$(hostname).mk

make all

# Check if the environment works
./bin/tp_testenv


# Run the code


for ((i = 10; i <= 100000; i *= 10));
do
	#./bin/tpPoisson1D_direct $i
	./bin/tpPoisson1D_iter 0 $i
	./bin/tpPoisson1D_iter 1 $i
	./bin/tpPoisson1D_iter 2 $i
done



# Clean
rm $(hostname).mk