#!/bin/bash

echo "Running test suite..."
commander test test/commander/hmmbuild_mult_seq_test.yaml
commander_1_exit=$?
commander test test/commander/dfam_tableizer_test.yaml
commander_2_exit=$?

exit_code=0

if [ $commander_1_exit -gt $commander_2_exit ]
then 
	exit_code=$commander_1_exit
else
	exit_code=$commander_2_exit
fi

echo "Finished"

exit $exit_code