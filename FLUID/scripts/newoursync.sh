#!/bin/bash -x

TO_TDS_FROM_BRIDGE=/home/fluid/TO_TDS_FROM_BRIDGE
TO_BRIDGE_FROM_TDS=/home/fluid/TO_BRIDGE_FROM_TDS
AGG_ISS_WORKSPACE_PATH=/home/ec2-user/data/AGG_ISS/WORKSPACE/workspace
COLAB_SHIM_WORKSPACE_PATH=/home/ec2-user/data/COLAB_SHIM/WORKSPACE/workspace
COLAB_SHIM=54.152.119.198
AGG_ISS=54.152.119.198
TDS=54.236.96.104
BRIDGE=52.5.52.204

#######################################
# run bridge on bridge
# run tds on tds
# run cnc on tds (for testing purposes)
#######################################

if [ $# -ne 2 ]
then
	echo "usage: $0 bridge|tds|cnc start|continue"
	exit 1
fi

MYHOSTNAME=$1
OPERATION=$2

if [ $MYHOSTNAME != "bridge" -a $MYHOSTNAME != "tds"  -a $MYHOSTNAME != "cnc" ]
then
	echo "usage: $0 bridge|tds|cnc"
	exit 1
fi



idempotentize() {
	d=$1
	if [ -d $d -a $OPERATION == 'start' ]
	then
		mv $d ~/DONE/$TS
		mkdir -p $d
	elif [ ! -d $d ] 
	then
		mkdir -p $d
	fi	
				
}




if [ $MYHOSTNAME == "bridge" ]
then
	TS=$(date +%s)
	mkdir -p ~/DONE/$TS

	idempotentize $TO_BRIDGE_FROM_TDS
	idempotentize $TO_TDS_FROM_BRIDGE

	while true
	do
		ssh fluid@${COLAB_SHIM} [ -d $COLAB_SHIM_WORKSPACE_PATH ] 
		colabshim_running=$(echo $?)
		ssh fluid@${TDS} [ -d $TO_TDS_FROM_BRIDGE -a -d $TO_BRIDGE_FROM_TDS ]
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
		# pull stuff from colab-shim
		rsync -a fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/request_path/ ${TO_TDS_FROM_BRIDGE}/
		rsync -a fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/response_path/ ${TO_TDS_FROM_BRIDGE}/
		# and push it to tds 
		rsync -a  ${TO_TDS_FROM_BRIDGE}/ fluid@${TDS}:${TO_TDS_FROM_BRIDGE}/

		# pull stuff from tds 
		rsync -a fluid@${TDS}:${TO_BRIDGE_FROM_TDS}/ ${TO_BRIDGE_FROM_TDS}/ 
		# and push it to colab-shim
		rsync -a ${TO_BRIDGE_FROM_TDS}/ fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/response_path
		rsync -a ${TO_BRIDGE_FROM_TDS}/ fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/request_path
		sleep 5 
	done


elif [ $MYHOSTNAME == "cnc" ]
then
	TS=$(date +%s)
	mkdir -p ~/DONE/$TS

	idempotentize $TO_TDS_FROM_BRIDGE
	idempotentize $TO_BRIDGE_FROM_TDS

	while true
	do
		# pull stuff from TDS 
		rsync -a fluid@${TDS}:${TO_TDS_FROM_BRIDGE}/ ${TO_TDS_FROM_BRIDGE}/
		# and push it to ISS
		rsync -a ${TO_TDS_FROM_BRIDGE}/ fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/response_path/
		rsync -a ${TO_TDS_FROM_BRIDGE}/ fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/request_path/

		# pull stuff from ISS 
		rsync -a fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/request_path/ ${TO_BRIDGE_FROM_TDS}/ 
		rsync -a fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/response_path/ ${TO_BRIDGE_FROM_TDS}/ 
		# and push it to TDS
		rsync -a ${TO_BRIDGE_FROM_TDS}/ fluid@${TDS}:${TO_BRIDGE_FROM_TDS}/
	done
		sleep 5

fi
