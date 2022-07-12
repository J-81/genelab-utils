def split_input_fastq_file(IDs_list, counts_dict, input_fastq, output_dir, prefix, suffix):

    # opening all potential output files (to avoid overhead of opening and closing on each entry)
    # dictionary to link cell IDs to their specific file handle
    files_dict = {}

    i = 0
    for cell_ID in IDs_list:

        i += 1

        curr_handle = "handle_" + str(i)

        files_dict[cell_ID] = curr_handle

        curr_file_path = str(output_dir) + str(prefix) + str(cell_ID) + str(suffix)

        files_dict[cell_ID] = open(curr_file_path, "w")


    unmatched_out = open(str(output_dir) + str(prefix) + "unmatched" + str(suffix), "w")

    # opening and running through our pooled input fastq
    with gzip.open(input_fastq, "rt") as fastq_in:

        # iterating fastq entries with biopython
        for header, seq, qual in FastqGeneralIterator(fastq_in):

            # looping through our cell IDs list so long as this entry hasn't been matched yet
            for cell_ID in IDs_list:

                # checking if ID is in header
                if "_" + str(cell_ID) + "_" in header:

                    files_dict[cell_ID].write("@%s\n%s\n+\n%s\n" % (header, seq, qual))

                    # adding a count for that cell ID to the counts dictionary
                    counts_dict[cell_ID] += 1

                    # stopping search for this entry
                    break

            else:

                # if not found, writing to unmatched output and counting
                unmatched_out.write("@%s\n%s\n+\n%s\n" % (header, seq, qual))

                counts_dict["unmatched"] += 1

    # closing all files
    for cell_ID in IDs_list:

        files_dict[cell_ID].close()

    unmatched_out.close()

    return(counts_dict)
