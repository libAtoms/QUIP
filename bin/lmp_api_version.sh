#!/bin/bash
#!/bin/bash
# Determine the LAMMPS API version from the git reported last changed date
# of file. The version is the UNIX timestap of the most recently changed file.

API_FILES="src/Potentials/quip_lammps_wrapper.f95"

QUIP_ROOT=$(dirname "$0")/..

# By default only show the date
VERBOSE=false

# Options for interactive use
while getopts 'v' opt; do
  case $opt in
    v)
      VERBOSE=true
      ;;
    \?)
      echo "Usage: $0"
      echo "'-v' to output name of newest file to stderr"
      exit 1
      ;;
  esac
done

function git_date
{
  if [ -e "$1" ]
  then
    DATE=$(git log -n1 --format="%at" "$1")
    [[ $? -eq 0 && ! -z $DATE ]] && echo "$DATE" || echo 0
  else
    echo 0
  fi
}

if [ -s "${QUIP_ROOT}/API_VERSION" ]; then
  echo -ne "$(cat "${QUIP_ROOT}/API_VERSION")"
  exit 0
elif [ -d "${QUIP_ROOT}/.git" ]; then
  # Sort everything as "DATE filename" list and take the last one
  API_VERSION=$(
    for I in $API_FILES
    do
      if [ -d "${QUIP_ROOT}/$(dirname "$I")" ]
      then
        cd "${QUIP_ROOT}/$(dirname "$I")"
        echo $(git_date $(basename "$I")) "$I"
        cd - >/dev/null
      fi
    done | sort -n | tail -1 )
  # filename to stderr if requested;
  # split string with parameter expansion
  if [ $VERBOSE == "true" ]; then
    echo "${API_VERSION##* }" >&2
  fi
  echo "${API_VERSION%% *}"
else
   echo "0"
   exit 0
fi
