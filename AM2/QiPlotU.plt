#set terminal postscript enhanced eps
set termopt enhanced
set autoscale
set title "Tip deflections of a composite box beam under sinesodial force"
set xlabel "TIME (sec)"
set ylabel "TIP DEFLECTIONS (in)"

plot "QiDisp_RK4.out" using 1:2 t 'u_1' with lines
