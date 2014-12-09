#!/bin/bash

# Usage:
# install: setup.sh -p packages 
# uninstall: setup.sh -u -p packages 

op=install

while getopts  "up:" flag
do
	case $flag in
		u) op=uinstall;;
		p) package=$OPTARG;;
	esac
done

if [ -e setup.log/$package.manifest ]; then
	for line in `cat setup.log/$package.manifest`; do
		if [ -f $line ]; then
			rm $line;
			echo "del $line";
		fi
	done
	rm setup.log/$package.manifest
	rm setup.log/$package.info
fi

geo_sim_sdk=http://10.76.5.56:8080/svn/geo_sim_sdk/packages
if [ $op == "install" ]; then
	svn export --force $geo_sim_sdk/$package/include include
	svn export --force $geo_sim_sdk/$package/Linux Linux
	svn export --force $geo_sim_sdk/$package/Windows Windows
	if [ ! -e setup.log ]; then
		mkdir setup.log;
	fi
	svn ls -R $geo_sim_sdk/$package/ \
		| grep -v "^doc/" | grep -v "/$" > setup.log/$package.manifest
	svn info $geo_sim_sdk/$package/ > setup.log/$package.info
fi
