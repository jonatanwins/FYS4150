import os


def read_last_n_lines(file_path, num_lines=100, chunk_size=4096):
    with open(file_path, "rb") as file:  # b for binary to not load the entire file
        file.seek(0, 2)  # Go to end of file
        file_size = file.tell()
        buffer = bytearray()
        lines = []
        pointer = file_size

        while pointer > 0 and len(lines) < num_lines:
            # Determine how much to read in this iteration
            pointer = max(
                0, pointer - chunk_size
            )  # move back a chunk at a time, but not before the beginning, 0
            file.seek(pointer)
            buffer = file.read(min(chunk_size, file_size))

            # Process the buffer to extract lines
            lines_in_chunk = buffer.split(b"\n")
            if len(lines_in_chunk) > 1:
                lines[:0] = [line.decode()[::-1] for line in lines_in_chunk[:-1]]
                buffer = lines_in_chunk[-1]  # remaining part of the last line
            else:
                buffer = lines_in_chunk[0]

        # Handle remaining buffer if necessary
        if buffer:
            lines.insert(0, buffer.decode()[::-1])

    return lines[:num_lines][::-1]


# adds suffix to new filename
def save_lines_to_file(file_path, lines, num_lines):
    base, ext = os.path.splitext(file_path)  # ext is .txt
    new_file_path = (
        f"{base}last_{num_lines}_lines{ext}"  # looks weird, but files always end with _
    )

    with open(new_file_path, "w") as new_file:
        new_file.write("\n".join(lines))

    return new_file_path


file_path = "code_p3/data/no_int_f_0_w_v_0_.txt"
