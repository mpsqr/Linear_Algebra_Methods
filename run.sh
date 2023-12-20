#!/bin/bash


cd tp/TP_Poisson_C_students_2022



cp ambre.mk ./$(hostname).mk

make all

# Check if the environment works
./bin/tp_testenv


# Run the code

./bin/tpPoisson1D_direct


./bin/tpPoisson1D_iter
./bin/tpPoisson1D_iter 1
./bin/tpPoisson1D_iter 2

# Clean
rm $(hostname).mk