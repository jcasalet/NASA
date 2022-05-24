#!/bin/bash


while true
do
	cp ISS/AGG/workspace/proto_path/* EARTH/COLAB/SHIM/workspace/proto_path/
	sleep 3
	cp EARTH/COLAB/SHIM/workspace/proto_path/* ISS/AGG/workspace/proto_path/
	sleep 3
done
