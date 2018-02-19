set terminal png font Helvetica 16
set output 'Results/RMSE/RMSE:0.025487/plot.png'
set xlabel 'Real data'
set ylabel 'Predicted data'
set xrange [-0.50:] 
 set yrange [-0.50:]
set style line 1 lc rgb '#0000A0' pt 7 ps 1 
set style line 2 lt 1 lc rgb 000000 lw 2 
f(x) = x 
plot 'Results/RMSE/RMSE:0.025487/regression' ls 1 title 'RMSE: 0.025487' ,  f(x) title 'Ideal' with lines ls 2, 'Results/RMSE/RMSE:0.025487/regression' using 1:2:($0+2) notitle with dots ls 1   
 