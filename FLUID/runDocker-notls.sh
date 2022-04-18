#!/bin/bash -x

# read these from a file, hardcode for now
DATA_PATH=/Users/jcasalet/Desktop/FLUID/CRISP/data
SCRIPT_PATH=/Users/jcasalet/Desktop/FLUID/DOCKER
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
	agg)
		HOSTNAME=agg
		IP_ADDR=192.168.56.101
		echo about to run container with role = $ROLE
		docker run --net $MYNET -h ${HOSTNAME} --rm  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro --ip ${IP_ADDR}  --add-host agg:192.168.56.101  --add-host colab1:192.168.56.102  --add-host colab2:192.168.56.103  ${IMAGE_NAME}  /scripts/runCrisp-notls.sh -r $ROLE
		;;
	colab1)
		HOSTNAME=colab1
		IP_ADDR=192.168.56.102
		echo "about to run container with role = $ROLE"
		docker run --net $MYNET -h ${HOSTNAME} --rm  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro  --ip ${IP_ADDR}  --add-host agg:192.168.56.101  --add-host colab1:192.168.56.102  --add-host colab2:192.168.56.103  ${IMAGE_NAME}  /scripts/runCrisp-notls.sh -r $ROLE 
		;;
	colab2)
		HOSTNAME=colab2
		IP_ADDR=192.168.56.103
		docker run --net $MYNET -h ${HOSTNAME} --rm  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro  --ip ${IP_ADDR}  --add-host agg:192.168.56.101  --add-host colab1:192.168.56.102  --add-host colab2:192.168.56.103  ${IMAGE_NAME}  /scripts/runCrisp-notls.sh -r $ROLE
		;;
	root)
		docker run -it --rm --user 0:0 $IMAGE_NAME /bin/bash
		;;
	*)
		echo "wrong usage: $ROLE"
		exit 1
		;;
esac
