#!/bin/bash

echo "Starting main script..."

declare -A server_files

# Mapping servers to the number of files they should process
server_files[salz1.idies.jhu.edu]=4
server_files[salz2.idies.jhu.edu]=4
server_files[salz3.idies.jhu.edu]=4
server_files[salz4.idies.jhu.edu]=3
server_files[salz6.idies.jhu.edu]=2
server_files[salz9.idies.jhu.edu]=4

processed_files=0

for server in "${!server_files[@]}"
do
    echo "Processing for server: $server"
    # Calculate the starting and ending files for the current server.
    start_file=$((processed_files + 1))
    end_file=$((processed_files + server_files[$server]))

    # Loop through the appropriate range of files for this server.
    for file_number in $(seq $start_file $end_file)
    do
        file="MANE_refseq_protein_chunk${file_number}.fasta"
        echo "    Running BLAST for file: $file"
        # Create the SSH command for clarity and debugging
        ssh_cmd="bash /home/choh1/intron_pos_conserv/blastdb/blast_run.sh /home/choh1/intron_pos_conserv/blastdb/$file 4 > /home/choh1/intron_pos_conserv/blastdb/blast_files/${file}_log.txt 2>&1 &"
        echo "    Executing: ssh $server \"$ssh_cmd\""
        ssh $server "$ssh_cmd"
        sleep 2
        ssh_exit_status=$?
        echo "    SSH exit status: $ssh_exit_status"
        if [ $ssh_exit_status -ne 0 ]; then
            echo "    Error executing SSH command on $server"
        fi
    done

    # Update the count of processed files.
    processed_files=$end_file
done

echo "Main script finished!"
