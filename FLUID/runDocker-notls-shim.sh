#!/bin/bash -x

# read these from a file, hardcode for now
DATA_PATH=/home/fluid/data
SCRIPT_PATH=/home/fluid/scripts
IMAGE_NAME=fluid
AGG_PORT=8888
AGG_EARTH_IP=100.27.49.222
AGG_ISS_IP=3.239.85.160
COLAB_EARTH_IP=3.235.222.119
COLAB_ISS_IP=3.239.83.239
COLAB_SHIM_IP=3.238.96.18

create_network() {
	MYNET=$1
	network_exists=$(docker network list | grep $MYNET)
	if [ -z "$network_exists" ]
	then
		docker network create --subnet=$SUBNET $MYNET
	fi
}

process_args() {
	ARGS=$1
	for arg in $ARGS
	do
		case $arg in
			-r |--role)
			ROLE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
			;;

       			*)
			echo "wrong usage: $ARGS"
			exit 1
       			;;
		esac
	done	
	echo $ROLE
}

ROLE=$(process_args $@)
echo my role in runDocker: $ROLE
#create_network $MYNET

case $ROLE in
	agg-iss)
		HOSTNAME=agg-iss
		;;
	agg-earth)
		HOSTNAME=agg-earth
		;;
	colab-iss)
		HOSTNAME=colab-iss
		;;
	colab-shim)
		HOSTNAME=colab-shim
		;;
	colab-earth)
		HOSTNAME=colab-earth
		;;
	root)
		docker run -it --rm --user 0:0 $IMAGE_NAME /bin/bash
		;;
	*)
		echo "wrong usage: $ROLE"
		exit 1
		;;
esac

echo about to run container with role = $ROLE

docker run -p 8888:8888 -h ${HOSTNAME}  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro  --add-host agg-iss:$AGG_ISS_IP --add-host agg-earth:$AGG_EARTH_IP --add-host colab-iss:$COLAB_ISS_IP --add-host colab-shim:$COLAB_SHIM_IP --add-host colab-earth:$COLAB_EARTH_IP --restart always  $IMAGE_NAME /scripts/runCrisp-notls-shim.sh -r $ROLE 
