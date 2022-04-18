#!/bin/bash -x

# read these from a file, hardcode for now
IMAGE_NAME=crisp
MYNET=mynet
SUBNET=192.168.56.0/24
AGG_PORT=8888
STATE_DIR=/data/STATE
WORKSPACE_DIR=/data/WORKSPACE
AGGRUNNING_FILE=aggrunning
ROUNDS=2

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

update_plans() {
	AGG_PORT=$1

	# update plan.yaml
	sed -i "s/.*agg_addr.*/    agg_addr: agg/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*agg_port.*/    agg_port: $AGG_PORT/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*disable_client_auth.*/    disable_client_auth: true/" ~/crisp/fl_plan/plan.yaml
	sed -i 's/.*disable_tls.*/    disable_tls: true/' ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*collaborator_count.*/    collaborator_count: 2/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*best_state_path.*/    best_state_path: \/data\/WORKSPACE\/crisp_best_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*init_state_path.*/    init_state_path: \/data\/WORKSPACE\/crisp_init_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*last_state_path.*/    last_state_path: \/data\/WORKSPACE\/crisp_last_.pbuf/" ~/crisp/fl_plan/plan.yaml
	sed -i "s/.*rounds_to_train.*/    rounds_to_train: $ROUNDS/" ~/crisp/fl_plan/plan.yaml
	
	cat ~/crisp/fl_plan/plan.yaml
	
	# update data.yaml
	echo "colab1,/data/col_0" > ~/crisp/fl_plan/data.yaml
	echo "colab2,/data/col_1" >> ~/crisp/fl_plan/data.yaml
	cat ~/crisp/fl_plan/data.yaml
	
	# update cols.yaml
	echo "collaborators:" > ~/crisp/fl_plan/cols.yaml	
	echo "- colab1" >> ~/crisp/fl_plan/cols.yaml	
	echo "- colab2" >> ~/crisp/fl_plan/cols.yaml	
	cat ~/crisp/fl_plan/cols.yaml

	# run federation_live_demo_agg.sh
	cd ~/crisp
	#./federation_live_demo_agg.sh

	# run commands from federation_live_demo_agg.sh for non-tls
	TEMPLATE=torch_cnn_mnist
	FED_WORKSPACE=workspace   # This can be whatever unique directory name you want
	FQDN=$(hostname --all-fqdns | awk '{print $1}')
	# 2) Create FL workspace
	rm -rf ${FED_WORKSPACE}
	fx workspace create --prefix ${FED_WORKSPACE} --template ${TEMPLATE}
	chmod 777 ${FED_WORKSPACE}
	cd ${FED_WORKSPACE}

	FED_DIRECTORY=`pwd`     # Get the absolute directory path for the workspace
	touch requirements.txt  ## PD: absolute hack, to remove/clear the requirements in the workspace :)

	# 3) Initialize our FL plan
	cp -r ../fl_plan/* plan/   # This is where all the collab detils should be
	cp -r ../fl_src/* src/
	fx plan initialize -a ${FQDN}

	# 4) Create certificate authority for workspace
	#fx workspace certify

	# 5) Create aggregator certificate
	#fx aggregator generate-cert-request --fqdn ${FQDN}

	# 6) Sign aggregator certificate
	#fx aggregator certify --fqdn ${FQDN} --silent # Remove '--silent' if you run this manually

	#7) Export the workspace so that it can be imported to the col nodes
	fx workspace export
	
	# scp workspace.zip
	#scp workspace/workspace.zip fluid@colab1:crisp
	#scp workspace/workspace.zip fluid@colab2:crisp
	cp workspace.zip ${STATE_DIR} 
	
}	

generate_cert_requests() {
	ROLE=$1
	cd ~/crisp
	mkdir workspace
	cp ${STATE_DIR}/workspace.zip workspace
	cd workspace
	unzip workspace.zip
	cat plan/plan.yaml
	#if [ $ROLE == colab1 ]
	#then	
	#	fx collaborator generate-cert-request -n $ROLE -d /data/col_0
	#else
	#	fx collaborator generate-cert-request -n $ROLE -d /data/col_1
	#fi
		
	#cp col_${ROLE}_to_agg_cert_request.zip ${STATE_DIR} 
	
}

activate_conda() {
	VENV=$1
	echo "activating conda env"
	conda init bash
	. ~/.bash_profile > /dev/null
	. ~/.bashrc > /dev/null
	conda activate $VENV 
}

run_agg() {
	if [ -d ${STATE_DIR} ]
	then
		rm -rf ${STATE_DIR}
	fi
	mkdir ${STATE_DIR}
	rm -rf ${WORKSPACE_DIR}/*
	update_plans 8888	
	touch ${STATE_DIR}/update_plans.done
	#while [[ ! -f ${STATE_DIR}/col_colab1_to_agg_cert_request.zip ]] || [[ ! -f ${STATE_DIR}/col_colab2_to_agg_cert_request.zip ]]
	#do
	#	sleep 5
	#done
	#cp ${STATE_DIR}/col_colab1_to_agg_cert_request.zip ~/crisp/workspace
	#cp ${STATE_DIR}/col_colab2_to_agg_cert_request.zip ~/crisp/workspace
	cd ~/crisp/workspace

	#fx collaborator certify -s --request-pkg col_colab1_to_agg_cert_request.zip
	#fx collaborator certify -s --request-pkg col_colab2_to_agg_cert_request.zip
	
	#zip -r cert.zip cert
	#cp cert.zip ${STATE_DIR} 

	fx aggregator start
	while [[ ! -d /data/WORKSPACE ]] || [[ ! -f /data/WORKSPACE/crisp_best_.pbuf ]]
	do
		sleep 5
	done
	cd ~/crisp/workspace
	fx model save -m /data/WORKSPACE/crisp_best_.pbuf
	
}

run_colab() {
	ROLE=$1
	while [ ! -f ${STATE_DIR}/update_plans.done ]
	do
		sleep 5
	done
	generate_cert_requests $ROLE
	#while [ ! -f ${STATE_DIR}/cert.zip ]
	#do
		#sleep 5
	#done
	#cp ${STATE_DIR}/cert.zip ~/crisp/workspace
	#cd ~/crisp/workspace
	#mv cert orig_cert
	#unzip cert.zip
	#mkdir cert/client
	#cp orig_cert/client/* cert/client	
	fx collaborator start -n $ROLE -p plan/plan.yaml -d plan/data.yaml
}
	

main() {
	activate_conda venv_3.7
	ROLE=$(process_args $@)
	echo my role in runCrisp: $ROLE
	case "$ROLE" in
		agg)
			run_agg	
			;;
	
		colab1|colab2)
			run_colab $ROLE	
			;;
		*)
			echo "wrong usage: $ROLE"
			exit 1
			;;
	
	esac
}

main $@ 
