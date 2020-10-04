cat farm-snapshot.txt | awk '{print $2}' | sort | uniq -c | sort -k1,1nr head -5 > top_5_users.txt

cat farm-snapshot.txt | grep "pathpip" | wc -l

grep "PEND" farm-snapshot.txt | wc -l

grep "RUN" farm-snapshot.txt | wc -l