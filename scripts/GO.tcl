set nEvents 10
set timestr 2025-10-15_16-58
set firmwareDir BI_DCT_FW
set outputDir DCT_raw_files
set firmwarePath ../$firmwareDir
set outputPath ./output/$outputDir

# channel mask
#          . layer2. layer1. layer0      
set eta1_1 000000000000000000000000
set eta1_2 000000000000000000000000
set eta1_3 000000000000000000000000
set eta1_4 000000000000000000000000
set eta1_5 000000000000000000000000
set eta1_6 000000000000000000000000

#          . layer2. layer1. layer0
set eta2_1 000000000000000000000000
set eta2_2 000000000000000000000000
set eta2_3 000000000000000000000000
set eta2_4 000000000000000000000000
set eta2_5 000000000000000000000000
set eta2_6 000000000000000000000000

set mask_eta1 $eta1_6$eta1_5$eta1_4$eta1_3$eta1_2$eta1_1
set mask_eta2 $eta2_6$eta2_5$eta2_4$eta2_3$eta2_2$eta2_1


open_hw
connect_hw_server
open_hw_target
current_hw_device [get_hw_devices xc7a200t_0]
refresh_hw_device -update_hw_probes false [lindex [get_hw_devices xc7a200t_0] 0]
set_property PROBES.FILE $firmwarePath/top.ltx [get_hw_devices xc7a200t_0]
set_property FULL_PROBES.FILE $firmwarePath/top.ltx [get_hw_devices xc7a200t_0]
set_property PROGRAM.FILE $firmwarePath/top.bit [get_hw_devices xc7a200t_0]
program_hw_devices [get_hw_devices xc7a200t_0]
refresh_hw_device [lindex [get_hw_devices xc7a200t_0] 0]

display_hw_ila_data [ get_hw_ila_data hw_ila_data_2 -of_objects [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]]
set_property CONTROL.DATA_DEPTH 128 [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]
set_property CONTROL.TRIGGER_POSITION 8 [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]
set_property TRIGGER_COMPARE_VALUE eq1'bR [get_hw_probes SMA_in_buf  -of_objects [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]]

set_property OUTPUT_VALUE_RADIX BINARY [get_hw_probes channel_mask_left -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ctrl_word_vio"}]]
set_property OUTPUT_VALUE $mask_eta1 [get_hw_probes channel_mask_left -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0]	-filter	{CELL_NAME=~"ctrl_word_vio"}]]
commit_hw_vio [get_hw_probes channel_mask_left -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0]	-filter	{CELL_NAME=~"ctrl_word_vio"}]]

set_property OUTPUT_VALUE_RADIX BINARY [get_hw_probes channel_mask_right -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ctrl_word_vio"}]]
set_property OUTPUT_VALUE $mask_eta2 [get_hw_probes channel_mask_right -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0]  -filter {CELL_NAME=~"ctrl_word_vio"}]]
commit_hw_vio [get_hw_probes channel_mask_right -of_objects [get_hw_vios -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ctrl_word_vio"}]]

#set now [clock seconds]
#set timestr [clock format $now -format "%y-%m-%d_%H-%M"]
#exec mkdir $outputDir/$timestr

for {set i 0} {$i < $nEvents} {incr i} {
    run_hw_ila [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]
    wait_on_hw_ila [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]
    display_hw_ila_data [upload_hw_ila_data [get_hw_ilas -of_objects [get_hw_devices xc7a200t_0] -filter {CELL_NAME=~"ila_elinks_inst"}]]
    write_hw_ila_data -legacy_csv_file -force -quiet $outputPath/$timestr/tmp_file_$i hw_ila_data_2
    if {($i % 10)==0} {
    	puts "Events: $i / $nEvents"
    }
}
close_hw_manager

exit

