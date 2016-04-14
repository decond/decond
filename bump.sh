#!/bin/bash

usage="bump.sh <version>"

if [ $# == 0 ] ; then
echo $usage
  exit 1;
fi

version="$1"

pyversion="python/decond/_version.py"
sed -i "/__version__/ s/\".*\"/\"$version\"/" $pyversion
sed -n '/__version__/ p' $pyversion
git add $pyversion

ftversion="fortran/src/varpars.F90"
sed -i "/decond_version/ s/\".*\"/\"$version\"/" $ftversion
sed -n '/decond_version/ p' $ftversion
git add $ftversion
