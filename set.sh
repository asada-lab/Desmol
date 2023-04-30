#!/bin/bash
ChemTS_folder=~/ChemTS-master
TARGET_folder=./gen10000
CURRENT_folder=$pwd
#
# copy original_ChemTS files.
#
mkdir -p $TARGET_folder
cp $ChemTS_folder/mcts_logp_improved_version/{cycle_scores.txt,logP_values.txt,SA_scores.txt,add_node_type.py,load_model.py,make_smile.py} $TARGET_folder
cp -r $ChemTS_folder/{data,RNN-model} $TARGET_folder
#
# Patch
#
cp mod_patch $TARGET_folder
cd $TARGET_folder
patch -p0 < mod_patch
cd ..
# 
#
#
echo 'Example: Run jupyter notebook. (gen10000.ipynb and ML.ipynb) in gen10000 folder.'
