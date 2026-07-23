#!/bin/bash

# Define remote and local variables
REMOTE_USER="vamika"
REMOTE_HOST="paramsmriti.nabi.res.in"
REMOTE_PORT="4422"
REMOTE_PATH="/home/vamika/bana/susy_new/susy/2d_Q16/susy/N_2_6x6_rt_g_0.6_1_hmc_PB_-1"
#REMOTE_PATH="/home/vamika/bana/susy_new/susy/2d_Q16/susy/update.c
LOCAL_PATH="/home/bana/susy/2d_Q16/susy"

# Extract the subdirectory name
SUBDIR_NAME=$(basename "$REMOTE_PATH")

# Print information
echo "Starting transfer from $REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH to $LOCAL_PATH..."

# Step 1: Determine a unique name for the new transfer
NEW_SUBDIR="$SUBDIR_NAME"

# If the directory already exists, append "_new" or "_new1", "_new2", etc.
COUNT=1
while [ -d "$LOCAL_PATH/$NEW_SUBDIR" ]; do
    NEW_SUBDIR="${SUBDIR_NAME}_new${COUNT}"
    ((COUNT++))
done

# Step 2: Transfer the directory from remote to local machine with new name
echo "Transferring directory from remote to local as '$NEW_SUBDIR'..."
scp -P "$REMOTE_PORT" -r "$REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH" "$LOCAL_PATH/$NEW_SUBDIR"

# Step 3: Check if the transfer was successful
if [ $? -eq 0 ]; then
    echo "Transfer completed successfully. New directory: $LOCAL_PATH/$NEW_SUBDIR"
else
    echo "Transfer failed. Please check your connection or input details."
fi

