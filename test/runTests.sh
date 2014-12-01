#!sh

for i in $(find . -perm +111 -type f)
do
./$i
done
