SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
for i in `cat hostfile`
do
  echo $i
  ssh $i "docker load -i ${ABS_DIRECTORY}/deepvariant.tar" &
done
