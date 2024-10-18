import os


def read_last_n_lines(file_path, num_lines=100):
    with open(file_path, "r") as file:
        # Read all lines and only keep the last num_lines lines
        lines = file.readlines()[-num_lines:]

    return [line.rstrip() for line in lines]  # Strip any extra newlines


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
