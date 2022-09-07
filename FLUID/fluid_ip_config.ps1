cd $home/.aws
$instances_file = ".\instances.json"
$hosts_file = ".\hosts"
$network_file = ".\fluid-network.sh"
$pem_file = "<<< FULLY QUALIFIED PEM FILE HERE. e.g. C:\TMP\JUPYTERLAB.PEM >>>>"
$fluid_node_names = @("bridge", "agg-earth", "colab-earth", "agg-iss", "colab-iss", "colab-shim")
$global:default_password = ""
$global:exception_password = @{"bridge" = "fluid"}   #list any nodes in the fluid network that do not use the default password. Make this an empty hashtable if there are not exceptions... = @{} 

function node_pwd {
	param (
        $node
    )
	if ($exception_password.$node -eq $null) {
		return $default_password
	} else {
		return ($exception_password.$node)
	}
}

echo "FLUID IP CONFIGURATION UTILITY"

# Generate json file of all AWS instances 
if (Test-Path $instances_file) {rm $instances_file}
$AccessKeyId = aws configure get aws_access_key_id
echo "Using AWS access key: $AccessKeyId"
aws ec2 describe-instances --query "Reservations[].Instances[].{Name: Tags[?Key=='Name']|[0].Value, PublicIP: PublicIpAddress, State: State.Name}"  --output json > $instances_file

#Abort if unable to access AWS instances
if ((!(Test-Path $instances_file)) -or ((Get-Item $instances_file).length -eq 0)) {
	Echo "Unable to access instances list in AWS. Please check credential configuration and retry."
	Exit
}

# create hashtable dictionary of AWS nodes that are part of the fluid architecture
$aws_instances = Get-Content $instances_file | ConvertFrom-Json
$fluid_nodes = $aws_instances.GetEnumerator() | ?{ $fluid_node_names -contains $_.Name } # filter out nodes that are not part of the fluid network

# Check that all fluid nodes are assigned a well-formed public IP address
$bad_ip = $fluid_nodes | ?{!($_.PublicIP -match '^((25[0-5]|(2[0-4]|1\d|[1-9]|)\d)(\.(?!$)|$)){4}$')}
if ($bad_ip) {
	echo "ERROR: All instances for the FLUID AWS architecture must have an assigned a Public IP."
	echo "Please activate the following node(s) and retry: `n"
	$bad_ip
	Exit
}

#Abort if any of the required fluid nodes have not been created
$missing_nodes = $fluid_node_names | foreach {if (!($fluid_nodes.Name -contains $_)) {$_}}
if ($missing_nodes) {
	echo "The following required node(s) have not been created. Please create and retry."
	write $missing_nodes
	Exit
}

echo "Generating IP configuration files for fluid nodes:"
echo $fluid_nodes

<#
	Create host file in format...
	
	127.0.0.1 localhost
	<public IP address of bridge> bridge     
	<public IP address of agg-earth> agg-earth
	<public IP address of agg-iss> agg-iss
	<public IP address of colab-earth> colab-earth
	<public IP address of colab-shim> colab-shim
	<public IP address of colab-iss> colab-iss
#>

"127.0.0.1 localhost" | Out-File -FilePath $hosts_file
$fluid_nodes | foreach {write $($_.PublicIP + " " + $_.Name)} | Out-File -FilePath $hosts_file -Append

<#

	Create fluid-network.sh file in format...

	DATA_PATH=~/data
	SCRIPT_PATH=~/scripts
	IMAGE_NAME=fluid
	AGG_PORT=8888
	AGG_EARTH_IP=100.27.49.222
	AGG_ISS_IP=3.239.85.160
	COLAB_EARTH_IP=3.235.222.119
	COLAB_ISS_IP=3.239.83.239
	COLAB_SHIM_IP=3.238.96.18
#>
write "DATA_PATH=/home/fluid/data" "SCRIPT_PATH=/home/fluid/scripts" "IMAGE_NAME=fluid" "AGG_PORT=8888" | Out-File -FilePath $network_file
$fluid_nodes | foreach {write $($_.Name.ToUpper().replace('-', '_')+"_IP" +"="+$_.PublicIP)} | Out-File -FilePath $network_file -Append

echo "`n`nIP configuration files created"

# copy hosts and fluid-network.sh files to each of the nodes in the fluid network
$default_password = Read-Host "`nEnter password for fluid accounts"
echo "`nCopying files to each fluid node."

$fluid_nodes | foreach {
	if ($_.State -eq "running") {
		Write "`n" $($_.Name + ": Copying files ...")
		$remote_dest = $("ec2-user@"+$_.PublicIP+":/tmp/hosts")
		scp -o "StrictHostKeyChecking no" -i $pem_file $hosts_file $remote_dest  # -o option to automatically add to known host if necessary
		$remote_dest = $("ec2-user@"+$_.PublicIP+":/tmp/fluid-network.sh")
		scp -i $pem_file $network_file $remote_dest
		# for Windows workstations, use iconv and sed to do the equivilant of dos2unix to fix the hosts and fluid-network.sh file
		$remote_dest = $("ec2-user@"+$_.PublicIP)
		$remote_commands = "iconv -f UTF-16LE -t UTF-8 /tmp/fluid-network.sh | sed 's/[^A-Za-z0-9_.;=~/]//g' > /tmp/tmp.sh; mv /tmp/tmp.sh /tmp/fluid-network.sh"
		ssh -i $pem_file $remote_dest $remote_commands
		$remote_commands = "iconv -f UTF-16LE -t UTF-8 /tmp/hosts | sed 's/[^A-Za-z0-9_.;=~\ -]//g' > /tmp/tmp_hosts; mv /tmp/tmp_hosts /tmp/hosts"
		ssh -i $pem_file $remote_dest $remote_commands
		# move the two files to their proper location
		$remote_commands = $("sudo mv /tmp/hosts /etc -f; echo " + $(node_pwd($_.Name)) + " | su - fluid -c 'cp /tmp/fluid-network.sh ~/scripts'; rm /tmp/fluid-network.sh -f")
		ssh -i $pem_file $remote_dest $remote_commands
	} else {
		Write "`n" $($_.Name + ": >> WARNING " + $_.Name + " had a public IP but was not running. Unable to copy files to that node.")
	}
}

echo "`n`nFLUID IP CONFIGURATION COMPLETE"

$openSSH = Read-Host "Do you want to open SSH sessions for the FLUID nodes now? (y/n)"
if ($openSSH -match "^y$|^yes$") {
	$fluid_nodes | foreach {
		if ($_.State -eq "running") {
			$sshCommand = $("-command ssh -i '"+$pem_file+"' ec2-user@"+$_.PublicIP)
			start powershell "-NoExit", $sshCommand
		} else {
			Write "`n" $($_.Name + ": >> WARNING will not open ssh for " + $_.Name + " because it is not running.")
		}
	}
}







