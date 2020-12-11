#!/bin/bash
sed -i '' "s/.*using namespace streamtube::model_.*/	using namespace streamtube::model_$1;/" "streamtube_concentration.cpp"
make streamtube_concentration
mv "streamtube_concentration" "../../bin/streamtube_concentration_$1"

