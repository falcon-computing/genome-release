#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export PATH=$DIR/bin:$PATH
export PATH=$DIR/tools/bin:$PATH
export LD_LIBRARY_PATH=$DIR/lib:$LD_LIBRARY_PATH
export LM_LICENSE_FILE=$DIR/license.lic
