-------------------------------------------------------------
HV
-------------------------------------------------------------
1, connect NHQ226L to gelab-serial04.mpp.mpg.de com1
set up the configurations
--> "Line" --> "Line 1 Configuration"
The serial interface is set to 9600Bit/s, 8Bit/character, no parity, 1Stop-Bit

2, set HV
ssh gelab@gelab-pc08
cd /remote/ceph/group/gedet/data/lab/2018/2018-06-22_8bfc8699_lm_regeb_in_new_lab/config_file
screen -r regeb_hv

venv legend-base
julia
include("communicate_with_nhq_226l_hv.jl")
get_nhq_226_status()
get_nhq_226_channel_1_hv()

---> check HV polarity chn1 REGeB -4500, chn2 GCDX +3200

---> to check the signal on oscilloscope, target HV -4500
set_nhq_226_channel_1_hv(100)  
set_nhq_226_channel_1_hv(500)  
set_nhq_226_channel_1_hv(4500)  

get_nhq_226_channel_1_hv()

GCDX
set_nhq_226_channel_2_hv(100)
set_nhq_226_channel_2_hv(500)
set_nhq_226_channel_2_hv(3200)


-------------------------------------------------------------
3, set DAQ
4, take data
-------------------------------------------------------------
screen -S regeb_daq
cd /remote/ceph2/group/gedet/data/lab/2018/2018-06-22_8bfc8699_lm_regeb_in_new_lab/raw_data
venv legend-base
//daqcore-scala-fadc -i ../config_file/regebchn1_gcdxchn2_daq_start_stop_bkg_20201007_reference.scala 
daqcore-scala-fadc -i ../config_file/gcdxchn1_daq_start_stop_bkg_20220118_reference.scala 
runFor(60.minutes, 61.minutes)



-------------------------------------------------------------
5, online monitor
-------------------------------------------------------------
cd /remote/ceph2/group/gedet/data/lab/2018/2018-06-22_8bfc8699_lm_regeb_in_new_lab/online_regeb_gcdx_monitor
venv legend-base
./start_regeb_gcdx_online_monitor.exe 

