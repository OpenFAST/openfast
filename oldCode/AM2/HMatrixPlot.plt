############
 # Title: [E3Tipu1]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H11_AM2.eps';
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
 set autoscale x
 set ylabel "DISPLACEMENT u_1 (inch)" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:2 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipu2]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H12_AM2.eps';
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
 set autoscale x
 set ylabel "DISPLACEMENT u_2 (inch)" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:3 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipu3]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H13_AM2.eps';
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
 set autoscale x
 set ylabel "DISPLACEMENT u_3 (inch)" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:4 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipp1]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H21_AM2.eps';
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
 set autoscale x
 set ylabel "ROTATION PARAMETER p_1" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:5 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipp2]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H22_AM2.eps';
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
 set autoscale x
 set ylabel "ROTATION PARAMETER p_2" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:6 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipp3]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H23_AM2.eps';
 unset key
 set border;
 set grid xtics lt 0;
 set grid ytics lt 0;
 set lmargin 11;
 set bmargin 4;
 set multiplot;
 set style line  1 lt  1 lw 3.0 pt  7 ps 1.5;
 unset title
 set xlabel "TIME (s)" font "Times,28";
 set autoscale x
 set ylabel "ROTATION PARAMETER p_3" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:7 with lines ls 1
 set nomultiplot;
 set output;
 
  ############
 # Title: [E3Tipp1]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H31_AM2.eps';
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
 set autoscale x
 set ylabel "ROTATION PARAMETER p_1" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:8 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipp2]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H32_AM2.eps';
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
 set autoscale x
 set ylabel "ROTATION PARAMETER p_2" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:9 with lines ls 1
 set nomultiplot;
 set output;
 
 ############
 # Title: [E3Tipp3]
 ############
 reset
 set terminal postscript eps enhanced color solid;
 set output 'H33_AM2.eps';
 unset key
 set border;
 set grid xtics lt 0;
 set grid ytics lt 0;
 set lmargin 11;
 set bmargin 4;
 set multiplot;
 set style line  1 lt  1 lw 3.0 pt  7 ps 1.5;
 unset title
 set xlabel "TIME (s)" font "Times,28";
 set autoscale x
 set ylabel "ROTATION PARAMETER p_3" font "Times,28";
 set autoscale y;
 plot "QiH_AM2.out" using 1:10 with lines ls 1
 set nomultiplot;
 set output;