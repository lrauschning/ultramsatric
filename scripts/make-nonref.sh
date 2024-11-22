#/bin/sh
for x in ./Pfam-B_*.ref.fa; do sed -e "s/-//g" $x > $(basename $x .ref.fa).fa; done
