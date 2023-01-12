#! /bin/bash
echo -e "starting Jupyter Lab... please wait!"
jupyter lab --no-browser --ip "*" >& jupyter-session.out &
sleep 6


node_id=$(hostname)
port=$(sed -n 's/.*localhost:\([0-9]*\)\/lab.*/\1/p' jupyter-session.out)

echo "command to use on your local machine:"
echo -e "ssh -L ${port}:${node_id}:${port} $USER@farm.cse.ucdavis.edu"

echo -e "and then use your web browser to open http://localhost:${port}"
