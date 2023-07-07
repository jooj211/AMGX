import os
import pandas as pd

def main():
    os.chdir("../build/csv")

    # read the CSV file into a pandas DataFrame
    df = pd.read_csv('summary_1.csv')

    # write the DataFrame to an Excel file
    df.to_excel('../sheets/summary_1.xlsx', index=False)

    # read the CSV file into a pandas DataFrame
    df = pd.read_csv('summary_2.csv')

    # write the DataFrame to an Excel file
    df.to_excel('../sheets/summary_2.xlsx', index=False)

    # read the CSV file into a pandas DataFrame
    df = pd.read_csv('summary_3.csv')

    # write the DataFrame to an Excel file
    df.to_excel('../sheets/summary_3.xlsx', index=False)

if __name__ == '__main__':
    main()
