#!/bin/bash

# Run the Python script
python run_ARDT_CyTRACK.py &

# Wait for all background processes to complete
wait

echo "Both scripts have completed."