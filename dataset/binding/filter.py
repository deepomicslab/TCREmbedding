def csv_filter(input_path='combined_dataset.csv', 
                 output_path='combined_dataset_filtered.csv',
                 epitope_length=20,
                 cdr_length=30):
    '''
    Modify /dataset/merged_data/combined_dataset.csv to Luu et al data format.

    :param input_path: input csv file path
    :param output_path: output csv file path
    :param epitope_length: maxium of epitope length
    :param cdr_length: maxium of cdr length
    :return: True if successful
    '''

    import pandas as pd

    print('Start reading csv file...')

    df = pd.read_csv(input_path)
    origin_rows = df.shape[0]
    print('Current rows of data: ' + str(df.shape[0]))
    
    # Feel free to modify the following statements if you want to use your own data.
    print('Start filtering csv file...')
    df = df[(df['Epitope'].str.match('^[A-Z]{1,' + str(epitope_length) + '}$')) & # Retain epitope sequences of capital letters with length less than the maximum 
            (~df['Epitope'].str.contains('B|J|O|U|X|Z')) & # Drop epitope sequences containing B/J/O/U/X/Z
            (df['CDR3'].str.match('^[A-Z]{1,'+ str(cdr_length) +'}$')) & # Retain CDR3 sequences of capital letters with length less than the maximum 
            (~df['CDR3'].str.contains('B|J|O|U|X|Z'))] # Drop CDR3 sequences containing B/J/O/U/X/Z
    filtered_rows = df.shape[0]
    
    print('Current rows of data: ' + str(filtered_rows) + 
          ', origin rows of data:' + str(origin_rows) + ', ' + 
          str(origin_rows - filtered_rows) + ' rows are removed.')
    
    # Remove duplicable data...
    print('Removing duplicable data...')
    df.drop_duplicates()
    
    print('Current rows of data: ' + str(df.shape[0]) + 
          ', previous rows of data:' + str(filtered_rows) + 
          ', ' + str(filtered_rows - df.shape[0]) + ' rows are removed. Totally ' + 
          str(origin_rows - df.shape[0]) + ' rows are removed.')
    
    print('Saving to csv file...')
    df.to_csv(path_or_buf=output_path, index=False)

    return True

# Example
csv_filter(input_path='combined_dataset.csv', 
            output_path='combined_dataset_filtered.csv',
            epitope_length=20,
            cdr_length=30)