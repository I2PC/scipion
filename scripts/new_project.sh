
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters. Just pass the project name."
    exit 128
fi

PROJECT_NAME=$1
$PWD
scipion python $SCIPION_HOME/scripts/create_project.py $PROJECT_NAME None $PWD
scipion project $PROJECT_NAME
