#!/bin/bash


cd tp/TP_Poisson_C_students_2022



cp ambre.mk ./$(hostname).mk

make all

# Check if the environment works
./bin/tp_testenv


# Run the code

./bin/tpPoisson1D_direct

# Clean
rm $(hostname).mk