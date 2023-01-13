#!/bin/bash -x

#configure the paths and IP addresses of the FLUID network
source ./crisp-config.sh


get_ip_addrs() {

	export AGG_EARTH_IP=$(grep agg-earth /etc/hosts | awk '{print $1}')	
	export AGG_ISS_IP=$(grep agg-iss /etc/hosts | awk '{print $1}')	
	export COLAB_EARTH_IP=$(grep colab-earth /etc/hosts | awk '{print $1}')	
	export COLAB_ISS_IP=$(grep colab-iss /etc/hosts | awk '{print $1}')	
	export COLAB_SHIM_IP=$(grep colab-shim /etc/hosts | awk '{print $1}')	


}


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

get_ip_addrs

#docker run -p 8888:8888 -h ${HOSTNAME}  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro  --add-host agg-iss:$AGG_ISS_IP --add-host agg-earth:$AGG_EARTH_IP --add-host colab-iss:$COLAB_ISS_IP --add-host colab-shim:$COLAB_SHIM_IP --add-host colab-earth:$COLAB_EARTH_IP --restart always  $IMAGE_NAME /scripts/runCrisp-notls-shim.sh -r $ROLE 

docker run -p 8888:8888 -h ${HOSTNAME}  --user=fluid:fluid  -v ${DATA_PATH}:/data:rw  -v ${SCRIPT_PATH}:/scripts:ro  --add-host agg-iss:$AGG_ISS_IP --add-host agg-earth:$AGG_EARTH_IP --add-host colab-iss:$COLAB_ISS_IP --add-host colab-shim:$COLAB_SHIM_IP --add-host colab-earth:$COLAB_EARTH_IP --restart no  $IMAGE_NAME /scripts/runCrisp-notls-shim.sh -r $ROLE 
