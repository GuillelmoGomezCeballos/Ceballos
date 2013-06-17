#!/bin/sh
kinit -f;
echo $KRB5CCNAME|awk '{split($1,a,"FILE:");print"cp "a[2]" ~/.krb5/ticket;"}' > kkk
cat kkk;
chmod a+x kkk;
./kkk;
rm -f kkk;
ls -l ~/.krb5/ticket;
