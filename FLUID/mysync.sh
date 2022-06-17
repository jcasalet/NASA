#!/bin/bash

rm FROM_EARTH/*
rm FROM_ISS/*

while true
do
	ssh fluid@colab-shim [ -d /home/fluid/data/WORKSPACE/workspace/proto_path ]
	earth_proto_exists=$(echo $?)
	ssh fluid@agg-iss [ -d /home/fluid/data/WORKSPACE/workspace/proto_path ]
	iss_proto_exists=$(echo $?)
	if [ $earth_proto_exists == 0 -a $iss_proto_exists == 0 ]
	then
		break
	else
		echo waiting for iss and earth nodes
		sleep 3
	fi
done
	

while true
do
	scp fluid@agg-iss:/home/fluid/data/WORKSPACE/workspace/proto_path/* FROM_ISS
	scp fluid@colab-shim:/home/fluid/data/WORKSPACE/workspace/proto_path/* FROM_EARTH
	scp FROM_EARTH/* fluid@agg-iss:/home/fluid/data/WORKSPACE/workspace/proto_path
	scp FROM_ISS/* fluid@colab-shim:/home/fluid/data/WORKSPACE/workspace/proto_path
	sleep 5 
done
