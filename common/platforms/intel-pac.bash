module_profile="aocl/17.1.1-pac"
license_dst="local"

# fpga info
export SW_BIT=$FALCON_HOME/fpga/sw.aocx
export SMEM_BIT=
export PMM_BIT=$FALCON_HOME/fpga/pmm.aocx

# platform tests enabled by default
export do_sw_test=0
export do_smem_test=0
export do_pmm_test=0
export do_blaze_test=1
