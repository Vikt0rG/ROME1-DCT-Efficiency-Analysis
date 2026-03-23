# source process_csv.sh iladata.txt
cat $1 | grep -v Sample | grep -v Rad | awk -F"," '{print $1,$4,$5}' > tmp.strip