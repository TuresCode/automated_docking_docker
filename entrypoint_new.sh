#!/bin/bash

   if [[ -x "/init.sh" ]]; then
     /init.sh
   fi

   # Source conda
   . /opt/conda/etc/profile.d/conda.sh

   # Activate the environment
   conda activate my_rdkit

   # Run the Python app
   python /flask_server/app.py

   # Execute any additional command passed to the script
   exec "$@"