set terminal postscript eps enhanced color font 'Helvetica,16'
set output 'convergence.eps'

set logscale x
set logscale y

set xrange [0.3:0.01] reverse

set format "%g"
#set format "%g"

plot 'convergence.txt' using (2/$1):2 w lp lw 2 title "k=0", \
     'convergence.txt' using (2/$1):3 w lp lw 2 title "k=1", \
     '-' using 1:2 w lp lw 2 title "ref. 1", \
     '-' using 1:2 w lp lw 2 title "ref. 1.5"
     0.1 0.8
     0.02 0.16
     e
     0.1 0.02
     0.02 0.00179
     e
