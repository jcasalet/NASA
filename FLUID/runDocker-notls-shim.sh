#!/bin/bash -x

# read these from a file, hardcode for now
DATA_PATH=/Users/jcasalet/Desktop/FLUID/CRISP/data
SCRIPT_PATH=/Users/jcasalet/Desktop/FLUID/DOCKER/SHIM
IMAGE_NAME=crisp
MYNET=mynet
SUBNET=192.168.56.0/24
AGG_PORT=8888
DATE=$(date +%s) 
AGGRUNNING_FILE=/tmp/aggrunning-$DATE

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
create_network $MYNET

case $ROLE in
	agg-iss)
		HOSTNAME=agg-iss
		IP_ADDR=192.168.56.101
		;;
	agg-earth)
		HOSTNAME=agg-earth
		IP_ADDR=192.168.56.102
		;;
	colab-iss)
		HOSTNAME=colab-iss
		IP_ADDR=192.168.56.103
		;;
	colab-earth)
		HOSTNAME=colab-earth
		IP_ADDR=192.168.56.104
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

docker run --net $MYNET -h ${HOSTNAME} --rm  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro --ip ${IP_ADDR}  --add-host agg-iss:192.168.56.101 --add-host agg-earth:192.168.56.102  --add-host colab-iss:192.168.56.103  --add-host colab-earth:192.168.56.104  ${IMAGE_NAME}  /scripts/runCrisp-notls-shim.sh -r $ROLE
