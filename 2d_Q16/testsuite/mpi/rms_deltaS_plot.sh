#!/bin/bash

parent_dir="/home/bana/susy/2d_Q16/susy/N_2_4x4_rt_g_4_0.5_hmc"
output_file="rms_deltaS_vs_eps.txt"
plot_file="rms_deltaS_loglog_eps.png"

rm -f "$output_file" "$plot_file"
echo "# epsilon(=1/nstep)   rms_deltaS" > "$output_file"

# ===============================
# Compute RMS(ΔS) after 100 configs
# ===============================
for folder in "$parent_dir"/nstep_*; do
    nstep=$(basename "$folder" | sed 's/nstep_//')
    input_file="$folder/data/deltaS.csv"

    if [[ -f "$input_file" ]]; then
        rms=$(awk -F'[,\t ]+' '
            NR > 101 && $2 ~ /^-?[0-9.]+$/ {
                sum += $2*$2
                count++
            }
            END {
                if (count > 0)
                    printf "%.6g", sqrt(sum/count)
                else
                    print "NaN"
            }' "$input_file")

        epsilon=$(awk -v n="$nstep" 'BEGIN { printf "%.6g", 1.0/n }')

        printf "%-15s %s\n" "$epsilon" "$rms" >> "$output_file"
    fi
done

echo "Saved RMS data → $output_file"

# ===============================
# Gnuplot with slope-2 reference
# ===============================
gnuplot << EOF
set terminal png size 800,600
set output "$plot_file"

set title "RMS(ΔS) vs ε"
set xlabel "ε"
set ylabel "RMS(ΔS)"

# Log10 axes
set logscale x 10
set logscale y 10

# Force numeric labels on major ticks
set format x "%g"
set format y "%g"

# Grid for both major and minor ticks
set grid xtics ytics mxtics mytics lw 1 lc rgb "#cccccc"

# Make sure tics are shown on x-axis
set xtics nomirror
set mxtics 10   # minor tics for logscale

set key top left

# Slope-2 reference line
stats "$output_file" using 1:2 nooutput
C = STATS_min_y / (STATS_min_x**2)
f(x) = C * x**2

plot \
  "$output_file" using 1:2 with linespoints pt 7 lw 2 title "RMS ΔS", \
  f(x) with lines lw 2 dt 2 title "∝ ε²"

EOF



echo "Plot saved as → $plot_file"

