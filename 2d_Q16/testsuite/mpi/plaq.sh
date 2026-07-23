#!/bin/bash

# Define the parent directory containing the folders
parent_dir="/home/bana/susy/2d_Q16/susy"

# Define the output files
output_file="avg_plaq.txt"
plot_file="avg_plaq_plot.png"

# Remove old files if they exist
[ -f "$output_file" ] && rm "$output_file"
[ -f "$plot_file" ] && rm "$plot_file"

# Initialize the output file
echo "# Col1 (6/rt_value)^2  Col2 (Average)" > "$output_file"

# Loop through all directories matching the pattern
for folder in "$parent_dir"/N_4_rt_g_*; do
    # Extract the value after "g_"
    rt_value=$(basename "$folder" | grep -oP '(?<=g_)[^_]+')

    # Calculate (rt_value / 6)^2
    col1_value=$(awk -v g="$rt_value" 'BEGIN {print (g / 6)^2}')

    # Define the input CSV file
    input_file="$folder/data/plaq.csv"

    # Check if the file exists
    if [[ -f "$input_file" ]]; then
        # Calculate the average of the sum of column 2 and column 3, starting from 10th row, skipping by 4
        avg=$(awk -F ',' 'NR >= 10 && (NR-10) % 5 == 0 {sum += ($2 + $3); count++} 
                          END {if (count > 0) print sum / count; else print "NaN"}' "$input_file")
        
        # Calculate reciprocal of col1_value
        reciprocal_col1=$(awk -v col1="$col1_value" 'BEGIN {print 1 / col1}')
        
        # Append the result to the output file
        printf "%-20s %s\n" "$reciprocal_col1" "$avg" >> "$output_file"

    else
        echo "File not found: $input_file" >&2
    fi
done

# Notify user
echo "Results saved to $output_file"

# Create a Gnuplot script for plotting
gnuplot_script="plot_script.gnuplot"

cat > "$gnuplot_script" <<EOF
set terminal png size 800,600
set output '$plot_file'
set title '<ssplaq + stplaq> vs 1/λ for U(4), 6x6 lattice'
set xlabel '1/λ'
set ylabel '<ssplaq + stplaq>'
set grid
plot '$output_file' using 1:2 with linespoints lc rgb 'blue' pt 7 title 'Data Points'
EOF

# Run Gnuplot to generate the plot
gnuplot "$gnuplot_script"

# Remove the Gnuplot script after plotting
rm "$gnuplot_script"

# Notify user
echo "Plot saved as $plot_file"

