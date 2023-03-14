#!/bin/bash -x

TO_TDS_FROM_BRIDGE=/home/fluid/TO_TDS_FROM_BRIDGE
TO_BRIDGE_FROM_TDS=/home/fluid/TO_BRIDGE_FROM_TDS
TO_ISS_FROM_TDS=/data/work/fluid/project/release/TO_ISS
TO_TDS_FROM_ISS=/data/work/fluid/project/release/TO_TDS
AGG_ISS_WORKSPACE_PATH=/home/ec2-user/data/AGG_ISS/WORKSPACE/workspace
COLAB_SHIM_WORKSPACE_PATH=/home/ec2-user/data/COLAB_SHIM/WORKSPACE/workspace
COLAB_SHIM=3.89.108.98
#TDS=192.48.188.134
TDS=44.207.6.238
BRIDGE=52.5.52.204
AGG_ISS=3.89.108.98
SLEEP_INTERVAL=5

#######################################
# run bridge on bridge
# run cnc on cnc (for testing purposes)
#######################################

if [ $# -ne 2 ]
then
	echo "usage: $0 bridge|cnc start|continue"
	exit 1
fi

MYHOSTNAME=$1
OPERATION=$2

if [ $MYHOSTNAME != "bridge" -a $MYHOSTNAME != "cnc" ]
then
	echo "usage: $0 bridge|cnc start|continue"
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
		ssh fluid@${TDS} [ -d $TO_ISS_FROM_TDS -a -d $TO_TDS_FROM_ISS ]
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
		# pull requests from tds 
		rsync -a fluid@${TDS}:${TO_TDS_FROM_ISS}/ ${TO_BRIDGE_FROM_TDS}/ 
		# and push requests to colab-shim
		rsync -a ${TO_BRIDGE_FROM_TDS}/ fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/request_path
		# pull responses from colab-shim
		rsync -a fluid@${COLAB_SHIM}:${COLAB_SHIM_WORKSPACE_PATH}/response_path/ ${TO_TDS_FROM_BRIDGE}/
		# and push responses to tds 
		rsync -a  ${TO_TDS_FROM_BRIDGE}/ fluid@${TDS}:${TO_ISS_FROM_TDS}/

		sleep $SLEEP_INTERVAL 
	done


elif [ $MYHOSTNAME == "cnc" ]
then
	TS=$(date +%s)
	mkdir -p ~/DONE/$TS

	idempotentize $TO_TDS_FROM_BRIDGE
	idempotentize $TO_BRIDGE_FROM_TDS

	while true
	do
		# pull responses from TDS 
		rsync -a fluid@${TDS}:${TO_ISS_FROM_TDS}/ ${TO_TDS_FROM_BRIDGE}/
		# and push responses to ISS
		rsync -a ${TO_TDS_FROM_BRIDGE}/ fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/response_path

		# pull requests from ISS 
		rsync -a fluid@${AGG_ISS}:${AGG_ISS_WORKSPACE_PATH}/request_path/ ${TO_BRIDGE_FROM_TDS}/ 
		# and push requests to TDS
		rsync -a ${TO_BRIDGE_FROM_TDS}/ fluid@${TDS}:${TO_TDS_FROM_ISS}/
		sleep $SLEEP_INTERVAL
	done

fi
