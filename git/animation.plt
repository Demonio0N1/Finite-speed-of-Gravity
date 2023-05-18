stats 'dat.dat' nooutput
set xrange [-50:50]
set yrange [-50:50]
set zrange [-50:50]


do for [i=1:int(STATS_blocks)] {
  splot for[j=1:2] "dat.dat" index (i-1) using 1:2:3 pt 7 lc rgb "black" notitle
  pause 0.0000001
}
