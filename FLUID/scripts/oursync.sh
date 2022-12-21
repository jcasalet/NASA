#!/bin/bash

TO_BRIDGE_FROM_COLABSHIM=/home/fluid/TO_BRIDGE_FROM_COLABSHIM
TO_COLABSHIM_FROM_BRIDGE=/home/fluid/TO_COLABSHIM_FROM_BRIDGE

TO_TDS_FROM_BRIDGE=/home/fluid/TO_TDS_FROM_BRIDGE
TO_BRIDGE_FROM_TDS=/home/fluid/TO_BRIDGE_FROM_TDS

TO_TDS_FROM_ISS=/home/fluid/TO_TDS_FROM_ISS
TO_ISS_FROM_TDS=/home/fluid/TO_ISS_FROM_TDS


if [ $# -ne 1 ]
then
	echo "usage: $0 bridge|tds"
	exit 1
fi

MYHOSTNAME=$1

if [ $MYHOSTNAME != "bridge" -a $MYHOSTNAME != "tds" ]
then
	echo "usage: $0 bridge|tds"
	exit 1
fi

if [ $MYHOSTNAME == "bridge" ]
then
	TS=$(date +%s)
	mkdir -p ~/DONE/$TS

	if [ -d $TO_BRIDGE_FROM_COLABSHIM ]
	then
		mv $TO_BRIDGE_FROM_COLABSHIM ~/DONE/$TS
	fi
	if [ -d $TO_COLABSHIM_FROM_BRIDGE ]
	then
		mv $TO_COLABSHIM_FROM_BRIDGE ~/DONE/$TS
	fi

	if [ -d $TO_BRIDGE_FROM_TDS ]
	then
		mv $TO_BRIDGE_FROM_TDS ~/DONE/$TS
	fi
	if [ -d $TO_TDS_FROM_BRIDGE ]
	then
		mv $TO_TDS_FROM_BRIDGE ~/DONE/$TS
	fi

	if [ -d $TO_TDS_FROM_ISS ]
	then
		mv $TO_TDS_FROM_ISS ~/DONE/$TS
	fi
	if [ -d $TO_ISS_FROM_TDS ]
	then
		mv $TO_ISS_FROM_TDS ~/DONE/$TS
	fi

	mkdir -p $TO_BRIDGE_FROM_COLABSHIM
	mkdir -p $TO_COLABSHIM_FROM_BRIDGE

	mkdir -p $TO_BRIDGE_FROM_TDS
	mkdir -p $TO_TDS_FROM_BRIDGE

	mkdir -p $TO_TDS_FROM_ISS
	mkdir -p $TO_ISS_FROM_TDS

	while true
	do
		ssh fluid@colab-shim [ -d /home/fluid/data/WORKSPACE/workspace/proto_path ]
		colabshim_running=$(echo $?)
		ssh fluid@tds [ -d $TO_TDS_FROM_ISS ]
		tds_running=$(echo $?)
		if [ $colabshim_running == 0 -a $tds_running == 0 ]
		then
			break
		else
			echo waiting for colab-shim and tds nodes
			sleep 3
		fi
	done
	

	while true
	do
		# get stuff from colab-shim
		rsync -a fluid@colab-shim:/home/fluid/data/WORKSPACE/workspace/proto_path/ $TO_BRIDGE_FROM_COLABSHIM
		# then push that to tds
		rsync -a $TO_BRIDGE_FROM_COLABSHIM/ fluid@tds:$TO_TDS_FROM_BRIDGE
		# get stuff from tds and push to colab-shim
		rsync -a $TO_BRIDGE_FROM_TDS/ fluid@colab-shim:/home/fluid/data/WORKSPACE/workspace/proto_path
		sleep 5 
	done

fi

if [ $MYHOSTNAME == "tds" ]
then
	TS=$(date +%s)
	mkdir -p ~/DONE/$TS

	if [ -d $TO_TDS_FROM_BRIDGE ]
	then
		mv $TO_TDS_FROM_BRIDGE ~/DONE/$TS
	fi

	if [ -d $TO_TDS_FROM_ISS ]
	then
		mv $TO_TDS_FROM_ISS ~/DONE/$TS
	fi

	mkdir -p $TO_TDS_FROM_BRIDGE
	mkdir -p $TO_TDS_FROM_ISS

	while true
	do
		ssh fluid@bridge [ -d $TO_BRIDGE_FROM_TDS ]
		bridge_running=$(echo $?)
		if [ $bridge_running == 0 ]
		then
			break
		else
			echo waiting for bridge node
			sleep 3
		fi
	done
	

	while true
	do
		# get stuff from bridge 
		rsync -a fluid@bridge:$TO_BRIDGE_FROM_COLABSHIM $TO_TDS_FROM_BRIDGE
		# and push it to ISS
		rsync -a $TO_TDS_FROM_BRIDGE/ fluid@agg-iss:/home/fluid/data/WORKSPACE/workspace/proto_path
		# get stuff from ISS
		rsync -a fluid@agg-iss:/home/fluid/data/WORKSPACE/workspace/proto_path/ $TO_TDS_FROM_ISS	
		rsync -a $TO_TDS_FROM_ISS/ fluid@bridge:$TO_BRIDGE_FROM_TDS
		sleep 5 
	done

fi
