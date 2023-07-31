cp basic_setup.sh /shared/

for i in `cat ../../hostfile`
do
  echo $i
  ssh $i "bash /shared/basic_setup.sh" &
done



##cp docker image to /shared/

#for i in `cat ../../hostfile`
# do
#  echo $i
#  ssh $i "docker load -i /shared/deepvariant.tar" &
#done
  
