using Sockets
using Base
using Formatting
#
# ##################################################
#
function get_nhq_226_hv()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   # following 3 lines remove the first polarity sign
   # no need, just skip the first letter
   # mytemp = collect(read_data) 
   # mytemp[1] = ' '
   # read_data=join(mytemp)  
   m = match(r"\+",read_data)
   plus_position =m.offset
   n = match(r"\-",read_data)
   minus_position =n.offset
   if  plus_position>0
      my_read_data=join([read_data[1:plus_position-1], "e",read_data[plus_position:end]])
   elseif  minus_position>0
      my_read_data=join([read_data[1:minus_position-1],"e",read_data[minus_position:end]])
   end
   close(nhqc)
   return my_read_data
end
#
# ##################################################
# 
function get_nhq_226_channel_1_hv()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   # following 3 lines remove the first polarity sign
   # no need, just skip the first letter
   # mytemp = collect(read_data)  
   # mytemp[1] = ' '
   # read_data=join(mytemp)
   m = match(r"\-",read_data)
   plus_position =m.offset
   n = match(r"\-",read_data)
   minus_position =n.offset
   if  plus_position>0
      my_read_data=join([read_data[1:plus_position-1], "e",read_data[plus_position:end]])
   elseif  minus_position>0
      my_read_data=join([read_data[1:minus_position-1],"e",read_data[minus_position:end]])
   end
   close(nhqc)
   return my_read_data
end
#
# ##################################################
# 
function get_nhq_226_channel_2_hv()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"U2\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   # following 3 lines remove the first polarity sign
   # no need, just skip the first letter
   # mytemp = collect(read_data)  
   # mytemp[1] = ' '
   # read_data=join(mytemp)
   m = match(r"\+",read_data)
   plus_position =m.offset
   n = match(r"\-",read_data)
   minus_position =n.offset
   if  plus_position>0
      my_read_data=join([read_data[1:plus_position-1], "e",read_data[plus_position:end]])
   elseif  minus_position>0
      my_read_data=join([read_data[1:minus_position-1],"e",read_data[minus_position:end]])
   end
   close(nhqc)
   return my_read_data
end
#
# ##################################################
#
function get_nhq_226_channel_3_hv()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10002)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   # following 3 lines remove the first polarity sign
   # no need, just skip the first letter
   # mytemp = collect(read_data)
   # mytemp[1] = ' '
   # read_data=join(mytemp)
   m = match(r"\+",read_data)
   plus_position =m.offset
   n = match(r"\-",read_data)
   minus_position =n.offset
   if  plus_position>0
      my_read_data=join([read_data[1:plus_position-1], "e",read_data[plus_position:end]])
   elseif  minus_position>0
      my_read_data=join([read_data[1:minus_position-1],"e",read_data[minus_position:end]])
   end
   close(nhqc)
   return my_read_data
end
#
# ##################################################
#
function get_nhq_226_channel_4_hv()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10002)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"U2\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   # following 3 lines remove the first polarity sign
   # no need, just skip the first letter
   # mytemp = collect(read_data)
   # mytemp[1] = ' '
   # read_data=join(mytemp)
   m = match(r"\+",read_data)
   plus_position =m.offset
   n = match(r"\-",read_data)
   minus_position =n.offset
   if  plus_position>0
      my_read_data=join([read_data[1:plus_position-1], "e",read_data[plus_position:end]])
   elseif  minus_position>0
      my_read_data=join([read_data[1:minus_position-1],"e",read_data[minus_position:end]])
   end
   close(nhqc)
   return my_read_data
end
#
# ##################################################
#
function set_nhq_226_hv(hvvalue::Real)
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   #---> read present HV value
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #println(read_data)
   #---> set HV ramp up speed
   sleep(0.5)
   write(nhqc,"V1=10\r\n")
   read_data = readline(nhqc)
   #---> set HV target value
   sleep(0.5)
   write(nhqc,"D1=$(round(hvvalue))\r\n") 
   read_data = readline(nhqc)
   #---> read present value
   sleep(0.5)
   write(nhqc,"G1\r\n") 
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   close(nhqc)
end
#
# ##################################################
#
function set_nhq_226_channel_1_hv(hvvalue::Real)
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   #---> read present HV value
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #println(read_data)
   #---> set HV ramp up speed
   sleep(0.5)
   write(nhqc,"V1=10\r\n")
   read_data = readline(nhqc)
   #---> set HV target value
   sleep(0.5)
   write(nhqc,"D1=$(round(hvvalue))\r\n")
   read_data = readline(nhqc)
   #---> read present value
   sleep(0.5)
   write(nhqc,"G1\r\n")
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   close(nhqc)
end
#
# ##################################################
#
function set_nhq_226_channel_2_hv(hvvalue::Real)
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   #---> read present HV value
   write(nhqc,"U2\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #println(read_data)
   #---> set HV ramp up speed
   sleep(0.5)
   write(nhqc,"V2=10\r\n")
   read_data = readline(nhqc)
   #---> set HV target value
   sleep(0.5)
   write(nhqc,"D2=$(round(hvvalue))\r\n")
   read_data = readline(nhqc)
   #---> read present value
   sleep(0.5)
   write(nhqc,"G2\r\n")
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   close(nhqc)
end
#
# ##################################################
#
function set_nhq_226_channel_3_hv(hvvalue::Real)
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10002)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   #---> read present HV value
   write(nhqc,"U1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #println(read_data)
   #---> set HV ramp up speed
   sleep(0.5)
   write(nhqc,"V1=10\r\n")
   read_data = readline(nhqc)
   #---> set HV target value
   sleep(0.5)
   write(nhqc,"D1=$(round(hvvalue))\r\n")
   read_data = readline(nhqc)
   #---> read present value
   sleep(0.5)
   write(nhqc,"G1\r\n")
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   close(nhqc)
end
#
# ##################################################
#
function set_nhq_226_channel_4_hv(hvvalue::Real)
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10002)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   #---> read present HV value
   write(nhqc,"U2\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #println(read_data)
   #---> set HV ramp up speed
   sleep(0.5)
   write(nhqc,"V2=10\r\n")
   read_data = readline(nhqc)
   #---> set HV target value
   sleep(0.5)
   write(nhqc,"D2=$(round(hvvalue))\r\n")
   read_data = readline(nhqc)
   #---> read present value
   sleep(0.5)
   write(nhqc,"G2\r\n")
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   close(nhqc)
end
#
#
# ##################################################
#
function get_nhq_226_status()
   scint_netcom_address = ip"134.107.13.65"
   nhqc = connect(scint_netcom_address, 10001)
   write(nhqc,"\r\n") # to assure synchronization
   readline(nhqc)
   readline(nhqc)
   write(nhqc,"S1\r\n")
   sleep(0.5)
   read_data = readline(nhqc)
   read_data = readline(nhqc)
   read_data = replace(read_data, "\r" => " ")
   read_data = replace(read_data, "\n" => " ")
   #  @printf("%s",read_data)
   close(nhqc)
   return read_data
end
