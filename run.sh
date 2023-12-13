#!/bin/bash


cd tp/TP_Poisson_C_students_2022



cp ambre.mk ./$(hostname).mk


make all

# Check if the environment works
./bin/tp_testenv


# Clean
rm $(hostname).mk