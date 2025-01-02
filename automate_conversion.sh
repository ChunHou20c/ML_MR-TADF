#!/usr/bin/env bash

# conversions
# this script is use to automate conversion process
# convert xyz file to mol

FILEPATH="./data/processed/"
for i in $(ls $FILEPATH)
do 
  basename $i .mol
  obabel $FILEPATH/$i -O "./data/smile/"$(basename $i .mol).smi 
done


# FILEPATH2="./data/molecules/raw/second_batch/"
# for i in $(ls $FILEPATH2)
# do 
#     basename $i .xyz
#     obabel $FILEPATH2/$i -O "./data/molecules/mol/second_batch/"$(basename $i .xyz).mol 
# done


# FILEPATH="./data/molecules/raw/first_batch/"
# for i in $(ls $FILEPATH)
# do 
#   basename $i .xyz
#   obabel $FILEPATH/$i -O "./data/molecules/smi/first_batch/"$(basename $i .xyz).smi 
# done
# 
# 
# FILEPATH2="./data/molecules/raw/second_batch/"
# for i in $(ls $FILEPATH2)
# do 
#     basename $i .xyz
#     obabel $FILEPATH2/$i -O "./data/molecules/smi/second_batch/"$(basename $i .xyz).smi
# done
