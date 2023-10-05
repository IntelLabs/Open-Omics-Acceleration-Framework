#!/bin/bash
# Input:
SECONDS=0
num_shards=$1


dir=$2

#cd $dir

candidate_files=`ls $dir/intermediate_results_dir*/candidates.tfrecord-*`

: '
for file in ${candidate_files}; do
    newname=$(echo "$file" | sed 's/-of-.....//'| sed 's/-//')

	echo "renaming $file to $newname"
	cp $file $newname

done
'


duration=$SECONDS
echo "Copy Part1: $duration sec"

SECONDS=0

# Total number of candidates across all shards
TotalCandidates=`awk '{ sum += $7 } END { print sum }' $dir/intermediate_results_dir*/*region.txt`

equal_distribution=`expr $TotalCandidates / $num_shards  `

echo "Total candidates across shards: $TotalCandidates equally distributing $equal_distribution candidates"

#echo ${equal_distribution} > $dir/max_candidates


duration=$SECONDS
echo "Total Part1: $duration sec"
SECONDS=0

regions_meta=`ls -v $dir/intermediate_results_dir*/*region.txt`
candidate_file_list=`ls -v $dir/intermediate_results_dir*/candidates.tfrecord-*.gz` #changed
list1=(${regions_meta})
list2=(${candidate_file_list})

# Check if the lists have the same length
: '
if [ ${#list1[@]} -ne ${#list2[@]} ]; then
  echo "Error: Lists have different lengths"
  exit 1
fi
'
echo "list1 ${#list1[@]}"
echo "list2 ${#list2[@]}"
duration=$SECONDS
echo "List: $duration sec"

SECONDS=0

rm -f $dir/master_list
# Iterate over the indices of one list (list1 in this example)
for ((i=0; i<${#list1[@]}; i++)); do
  item1="${list1[i]}"
  item2="${list2[i]}"
  #start=$(date +%s%3N) 
  
  #candidates_in_shards=`awk '{ sum += $7 } END { print sum }' $item1`
  
  sum=0
  while IFS= read -r line
  do
	  # Convert line into array
	  arr=($line)

	  # Get the seventh word
	  word=${arr[6]}

	  # Add to sum
	  sum=$((sum + word))
	done < $item1 #"yourfile.txt"
   #echo $sum
  candidates_in_shards=$sum   
  
  
  
  
  
  #end=$(date +%s%3N)
  #duration=$((end-start))
  #echo "Find: $duration millisec"
  #echo "List1 item: $item1 List2 item: $item2 No of Candidates: ${candidates_in_shards}"
  #SECONDS=0
  #start=$(date +%s%3N)
  echo "$item2 $item1 ${candidates_in_shards}" >> $dir/master_list 
  #end=$(date +%s%3N)
	
  #duration=$SECONDS
  #duration=$((end-start))
  #echo "Find1: $duration millisec"
done

duration=$SECONDS
echo "Part1: $duration sec"
SECONDS=0
rm -rf get_range
g++ get_range.cpp -o get_range

rm -f $dir/candidate_file_list $dir/counter_file_list $dir/offset_list $dir/max_candidates

./get_range $dir/master_list ${equal_distribution} $dir/candidate_file_list $dir/counter_file_list $dir/offset_list $dir/max_candidates

: '
ls $dir/candidate_file_list $dir/counter_file_list $dir/offset_list $dir/max_candidates

cat $dir/candidate_file_list 
cat $dir/counter_file_list 
cat $dir/offset_list
cat $dir/max_candidates
'
duration=$SECONDS
echo "Part2: $duration sec"

#3. offset list
#4. Region list
#4. max_candidates
