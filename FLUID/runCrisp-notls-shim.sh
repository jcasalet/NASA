#!/bin/bash -x

# read these from a file, hardcode for now
IMAGE_NAME=crisp
MYNET=mynet
SUBNET=192.168.56.0/24
STATE_DIR=/data/STATE
WORKSPACE_ISS_DIR=/data/WORKSPACE/ISS
WORKSPACE_EARTH_DIR=/data/WORKSPACE/EARTH
ROUNDS=2
AGG_ISS_HOST=agg-iss
AGG_ISS_PORT=8889
AGG_EARTH_HOST=agg-earth
AGG_EARTH_PORT=8888

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

update_plans_iss() {
	# update plan.yaml
	sed -i "s/.*agg_addr.*/    agg_addr: $AGG_ISS_HOST/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*agg_port.*/    agg_port: $AGG_ISS_PORT/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*disable_client_auth.*/    disable_client_auth: true/" ~/crisp/fl_plan/plan.yaml
	sed -i 's/.*disable_tls.*/    disable_tls: true/' ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*collaborator_count.*/    collaborator_count: 1/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*best_state_path.*/    best_state_path: ${WORKSPACE_ISS_DIR}\/crisp_best_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*init_state_path.*/    init_state_path: ${WORKSPACE_ISS_DIR}\/crisp_init_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*last_state_path.*/    last_state_path: ${WORKSPACE_ISS_DIR}\/crisp_last_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*rounds_to_train.*/    rounds_to_train: $ROUNDS/" ~/crisp/fl_plan/plan.yaml
	
	cat ~/crisp/fl_plan/plan.yaml
	
	echo "colab-iss,/data/col_1" > ~/crisp/fl_plan/data.yaml
	cat ~/crisp/fl_plan/data.yaml
	
	echo "collaborators:" > ~/crisp/fl_plan/cols.yaml	
	echo "- colab-iss" >> ~/crisp/fl_plan/cols.yaml	
	cat ~/crisp/fl_plan/cols.yaml

	cd ~/crisp
	rm -rf workspace 
	fx workspace create --prefix workspace --template torch_cnn_mnist 
	chmod 777 workspace
	cd workspace 

	cp -r ../fl_plan/* plan/
	cp -r ../fl_src/* src/
	fx plan initialize -a $AGG_ISS_HOST
}	

update_plans_earth() {
	# update plan.yaml
	sed -i "s/.*agg_addr.*/    agg_addr: $AGG_EARTH_HOST/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*agg_port.*/    agg_port: $AGG_EARTH_PORT/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*disable_client_auth.*/    disable_client_auth: true/" ~/crisp/fl_plan/plan.yaml
	sed -i 's/.*disable_tls.*/    disable_tls: true/' ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*collaborator_count.*/    collaborator_count: 1/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*best_state_path.*/    best_state_path: ${WORKSPACE_EARTH_DIR}\/crisp_best_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*init_state_path.*/    init_state_path: ${WORKSPACE_EARTH_DIR}\/crisp_init_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*last_state_path.*/    last_state_path: ${WORKSPACE_EARTH_DIR}\/crisp_last_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*rounds_to_train.*/    rounds_to_train: $ROUNDS/" ~/crisp/fl_plan/plan.yaml
	
	cat ~/crisp/fl_plan/plan.yaml
	
	echo "colab-earth,/data/col_0" > ~/crisp/fl_plan/data.yaml
	cat ~/crisp/fl_plan/data.yaml
	
	echo "collaborators:" > ~/crisp/fl_plan/cols.yaml	
	echo "- colab-earth" >> ~/crisp/fl_plan/cols.yaml	
	cat ~/crisp/fl_plan/cols.yaml

	cd ~/crisp
	rm -rf workspace 
	fx workspace create --prefix workspace --template torch_cnn_mnist 
	chmod 777 workspace
	cd workspace 

	cp -r ../fl_plan/* plan/
	cp -r ../fl_src/* src/
	fx plan initialize -a $AGG_EARTH_HOST 
}

activate_conda() {
	VENV=$1
	echo "activating conda env"
	conda init bash
	. ~/.bash_profile > /dev/null
	. ~/.bashrc > /dev/null
	conda activate $VENV 
}

run_agg_iss() {
	ROLE=agg-iss
	if [ -d $WORKSPACE_ISS_DIR ]
	then
		rm -rf $WORKSPACE_ISS_DIR
	fi
	mkdir $WORKSPACE_DIR
	update_plans_iss agg-iss 
	fx aggregator shim
}

run_agg_earth() {
	ROLE=agg-earth
	if [ -d $WORKSPACE_EARTH_DIR ]
	then
		rm -rf $WORKSPACE_EARTH_DIR
	fi
	mkdir $WORKSPACE_DIR
	update_plans_iss agg-iss 
	fx aggregator start
}

run_colab_iss() {
	ROLE=colab-iss
	update_plans_iss 
	fx collaborator start -n $ROLE -p plan/plan.yaml -d plan/data.yaml
}

run_colab_earth() {
	ROLE=colab-earth
	update_plans_earth 
	fx collaborator shim  -n $ROLE -p plan/plan.yaml -d plan/data.yaml
}
	

main() {
	activate_conda venv_3.7
	ROLE=$(process_args $@)
	HOST=$(echo $ROLE | cut -d- -f1)
	LOC=$(echo $ROLE | cut -d- -f2)
	echo my role in runCrisp: $ROLE
	case "$ROLE" in
		agg-iss)
			run_agg_iss
			;;
	
		agg-earth)
			run_agg_earth	
			;;
	
		colab-iss)
			run_colab_iss
			;;
		colab-earth)
			run_colab_earth
			;;
		*)
			echo "wrong usage: $ROLE"
			exit 1
			;;
	
	esac
}

main $@ 
