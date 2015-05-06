do for [i=3:11] {
set terminal png size 400,300 enhanced font "Helvetica,10"
unset key
set xlabel "<m>"
set ylabel "T"
outfile = sprintf('/home/prozek/NxN/mvsT'.i.'.png',i)
set output outfile
plot '/home/prozek/NxN/output'.i.'.csv' using 1:3 with lines
}
