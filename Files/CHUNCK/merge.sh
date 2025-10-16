count=0
for f in chunk_*; do
  root -b -q "MergeAllSparse.C(\"$f\",\"merged_$f.root\")" &
  ((count++))
  if ((count % 4 == 0)); then  # run 4 at a time
    wait
  fi
done
wait

