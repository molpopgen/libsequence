#!/bin/bash

#Process the header
#pandoc -S -s -c pandoc.css header.md -o header.html

 cd md

 for i in *.md
 do
     n=`basename $i .md`
     echo processing $i
     pandoc -S -s -c pandoc.css  -o ../$n.html $i
 done

cd ..

# cd doc

# for i in *.md
# do
#     n=`basename $i .md`
#     echo processing $i
#     pandoc -S -s -c pandoc.css  -o ../$n.html $i
# done

#cd ../blogs

#for i in *.md; 
#do
#	pandoc -S -s -c ../pandoc.css -H ../header.html -o ./${i//md/html} $i	
#	git add ./${i//md/html}
#done

#git commit -am "update"; 
#git push
