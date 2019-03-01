#!/bin/bash

export FALCON_DIR=$FALCON_HOME
export FCSBIN=$FALCON_DIR/bin/fcs-genome
export BWABIN=$FALCON_DIR/tools/bin/bwa-flow
export MMAPBIN=$FALCON_DIR/tools/bin/minimap-flow
export BLAZEBIN=$FALCON_DIR/blaze/bin/nam
export GATK3=$FALCON_DIR/tools/package/GATK3.jar
export GATK4=$FALCON_DIR/tools/package/GATK4.jar

# no validation in the env setup script
