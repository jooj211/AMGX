import os
import re
import sys
import csv
import subprocess

def parse_file(file_path):
    data = {}
    method = os.path.splitext(os.path.splitext(os.path.basename(file_path))[0])[0]
    data['method'] = method
    
    with open(file_path, 'r') as f:
        content = f.read()
        data['number_of_points'] = re.search('numeroPontos: (.*)\n', content).group(1).strip() if re.search('numeroPontos: (.*)\n', content) else None
        data['levels'] = re.search('Number of Levels: (.*)\n', content).group(1).strip() if re.search('Number of Levels: (.*)\n', content) else None
        data['total_memory_usage'] = re.search('Total Memory Usage: (.*) GB\n', content).group(1).strip() if re.search('Total Memory Usage: (.*) GB\n', content) else None
        data['total_iterations'] = re.search('Total Iterations: (.*)\n', content).group(1).strip() if re.search('Total Iterations: (.*)\n', content) else None
        data['final_residual'] = re.search('Final Residual: (.*)\n', content).group(1).strip() if re.search('Final Residual: (.*)\n', content) else None
        data['total_time'] = re.search('Total Time: (.*)\n', content).group(1).strip() if re.search('Total Time: (.*)\n', content) else None
       
    return data

def main():
    os.chdir("../build")
    output_path = './output_' + sys.argv[1]
    output_files = os.listdir(output_path)

    parsed_data = []

    for file in output_files:
        if file.endswith('.txt'):
            parsed_data.append(parse_file(os.path.join(output_path, file)))

    os.chdir("csv")

    with open('summary_' + sys.argv[1] + '.csv', 'w', newline='') as csvfile:
        fieldnames = ['method', 'number_of_points', 'levels', 'total_memory_usage', 'total_iterations', 'final_residual', 'total_time']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in parsed_data:
            writer.writerow(data)

if __name__ == '__main__':
    main()
    os.chdir("../../python")
    command = ['python', 'makeSheets.py', sys.argv[1]]
    subprocess.run(command)
