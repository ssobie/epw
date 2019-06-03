#!/bin/bash

scenario='rcp85' 
lon=-123.234678
lat=49.253745
infile='CAN_BC_UBCWesbrook_offset_from_POINT-ATKINSON_1106200_CWEC.epw'
name='UBCWesbrook'
qsub -N "epw.roll.DW" -v scenario=$scenario,method='roll',rlen=21,lon=$lon,lat=$lat,infile=$infile,name=$name run.epw.morph.pbs
qsub -N "epw.seas.DW" -v scenario=$scenario,method='seasonal',rlen='0',lon=$lon,lat=$lat,infile=$infile,name=$name run.epw.morph.pbs


