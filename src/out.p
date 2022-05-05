set terminal png size 1920,1080
set output 'graph.png'
set title "Bureau d'Etude"


set linetype 4 dt 1
set style line 1 lt 2 lc rgb "orange" lw 2


plot 'resultats/MyFile.csv' ls 1 