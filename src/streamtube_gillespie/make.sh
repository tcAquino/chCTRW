#!/bin/bash
sed -i '' "s/.*using namespace streamtube::model_.*/	using namespace streamtube::model_$1;/" "streamtube_gillespie.cpp"
make streamtube_gillespie
mv "streamtube_gillespie" "../../bin/streamtube_gillespie_$1"

