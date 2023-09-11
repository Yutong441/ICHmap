#! /bin/bash
cd /opt
mkdir -p data/.julia_lib
export PATH=$PATH:/opt/.julia/bin
export JULIA_DEPOT_PATH=/opt/data/.julia_lib
julia ICHmap/install_package.jl
