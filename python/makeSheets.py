import sys
import os
import pandas as pd

def main():
    os.chdir("../build/csv")

    # read the CSV file into a pandas DataFrame
    df = pd.read_csv('summary_' + sys.argv[1] + '.csv')

    # write the DataFrame to an Excel file
    df.to_excel('../sheets/summary_' + sys.argv[1] + '.xlsx', index=False)

if __name__ == '__main__':
    main()
