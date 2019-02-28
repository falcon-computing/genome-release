module_profile="sdx/17.4"
license_dst="hwc"

# fpga info
export SW_BIT=$FALCON_HOME/fpga/sw.xclbin
export SMEM_BIT=$FALCON_HOME/fpga/sw.xclbin
export PMM_BIT=$FALCON_HOME/fpga/pmm.xclbin

# platform tests enabled by default
export do_sw_test=1
export do_smem_test=1
export do_pmm_test=1
export do_blaze_test=1
