import os


def read_last_n_lines(file_path, num_lines=100):
    with open(file_path, "r") as file:
        # Read all lines and only keep the last num_lines lines
        lines = file.readlines()[-num_lines:]

    return [line.rstrip() for line in lines]  # Strip any extra newlines


# adds suffix to new filename
def save_lines_to_file(file_path, lines, num_lines, target_folder):
    base, ext = os.path.splitext(os.path.basename(file_path))  # ext is .txt
    new_file_path = os.path.join(
        target_folder, f"{base}last_{num_lines}_lines{ext}"
    )  # looks weird but all files end in _

    with open(new_file_path, "w") as new_file:
        new_file.write("\n".join(lines))

    return new_file_path


def process_folder(source_folder, target_folder, num_lines=100):
    # make target folder unless it exists
    os.makedirs(target_folder, exist_ok=True)

    for file_name in os.listdir(source_folder):
        file_path = os.path.join(source_folder, file_name)

        if os.path.isfile(file_path):
            lines = read_last_n_lines(file_path, num_lines)
            new_file_path = save_lines_to_file(
                file_path, lines, num_lines, target_folder
            )
            print(f"Last {num_lines} lines saved to: {new_file_path}")


source_folder = "data"
num_lines = 100
target_folder = f"last_{num_lines}_lines_w_v_from_2.15_to_2.25_step_0.01_no_int"

process_folder(source_folder, target_folder, num_lines)
