import os


def read_last_n_lines(file_path, num_lines=100):
    with open(file_path, "rb") as file:  # b for binary to not load the entire file
        file.seek(0, 2)  # end of file
        file_size = file.tell()  # tell is number of bytes to postion
        buffer = bytearray()  # mutable
        lines = []
        pointer = file_size - 1  # start of last byte

        # Read backwards and collect lines
        while pointer >= 0 and len(lines) < num_lines:
            file.seek(pointer)
            buffer.extend(file.read(1))

            # If we hit a newline
            if buffer[-1:] == b"\n":
                line = buffer[:-1].decode()[::-1]
                lines.append(line)
                buffer.clear()

            pointer -= 1

        # if no newline, append at the beginning
        if buffer:
            lines.append(buffer.decode()[::-1])

    return lines[::-1]


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
num_lines = 100
lines = read_last_n_lines(file_path, num_lines)
last_100_lines_file = save_lines_to_file(file_path, lines, num_lines)
print(f"Last {num_lines} lines saved to: {last_100_lines_file}")
