import os
import subprocess

executable = "examples/meu_teste2"
method_flag = "-c"
method_folder = "../src/configs"
step_flag = "-s"
step_size = "0.0002"
output_folder = "output_3"

method_files = os.listdir(method_folder)
total_files = len(method_files)

# Change the working directory to "../build"
os.chdir("../build")

# Create the output directory if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Counter to track progress
counter = 0

for method_file in method_files:
    method_path = os.path.join(method_folder, method_file)
    command = [executable, method_flag, method_path, step_flag, step_size]

    # Define the output file path within the output folder
    output_file = os.path.join(output_folder, f"{method_file}.txt")

    # Execute the command and redirect the output to the file
    with open(output_file, "w") as file:
        subprocess.run(command, stdout=file, stderr=subprocess.STDOUT)

    # Increment the counter and calculate the progress percentage
    counter += 1
    progress = (counter / total_files) * 100
    print(f"Progress: {progress:.2f}%")
