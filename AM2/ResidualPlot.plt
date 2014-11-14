 reset
 set terminal postscript eps enhanced color solid;
 set output 'Residual_TS1.eps';
 unset key
 set border;
 set grid xtics lt 0;
 set grid ytics lt 0;
 set lmargin 10;
 set bmargin 4;
 set multiplot;
 set style line  1 lt  1 lw 3.0 pt  7 ps 1.5;
 unset title
 set xlabel "TIME (s)" font "Times,28";
 set logscale;
 set autoscale x
 set ylabel "DISPLACEMENT u_1 (inch)" font "Times,28";
 set autoscale y;
 plot "Residual_TimeStep1.txt" using 1:2 with lines ls 1
 set nomultiplot;
 set output;