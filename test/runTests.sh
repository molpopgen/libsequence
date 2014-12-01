#!sh

for i in $(find . -executable -type f)
do
./$i
done
