#! /bin/bash
echo -e "starting Jupyter Lab... please wait!\n"
jupyter-lab --no-browser --ip "*" --notebook-dir ~/ >& jupyter-session.out &
sleep 15


node_id=$(hostname)
port=$(sed -n 's/.*localhost:\([0-9]*\)\/lab.*/\1/p' jupyter-session.out)

echo -e "1) copy-paste the following command in your local terminal\n"

GREEN='\033[0;32m'
NOCOLOR='\033[0m'
echo -e "${GREEN}ssh -L ${port}:${node_id}:${port} $USER@farm.cse.ucdavis.edu${NOCOLOR}\n"
echo -e "2) In your internet browser, open the following URL ${GREEN}http://localhost:${port}"
