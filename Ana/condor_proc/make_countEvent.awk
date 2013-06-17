BEGIN {
}
{
if(NR!=1) printf(",\n");
printf("   \"old/histo_%s_all_noskim.root\"",$5);
}
END{
printf("\n");
}
