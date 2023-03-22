import pandas as pd

def lista_para_string(lst):
    return ','.join(map(str, lst))

def tag_elements(tax_file, log):

    df = pd.read_csv(tax_file, sep='\t')

    # Tratamento para adicionar colunas para aplicação da tag
    df['Element-ID'] = df['Element-ID'].replace('.*/','', regex=True)
    df['start'] = df['Element-ID'].replace('-.*','', regex=True)
    df['start'] = df['start'].replace('.*:', '', regex=True)
    df['end'] = df['Element-ID'].replace('.*-', '', regex=True)
    df['contig'] = df['Element-ID'].replace(':.*', '', regex=True)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    # Indentificando linhas que existe overlap
    matched_elements = []
    for i, row in df.iterrows():
        matched_ids = []
        filtered_df = df[df['contig'] == row['contig']]
        for j, other_row in filtered_df.iterrows():
            if i != j and other_row['start'] <= row['end'] + 50 and other_row['end'] >= row['start'] - 50 and other_row['Family'] != row['Family']:
                matched_ids.append(other_row['Element-ID'])
        matched_elements.append(matched_ids)

    # Atribuir a lista matched_elements à coluna 'Overlaped_Element_ID'
    df['Overlaped_Element_ID'] = matched_elements

    # Tratamento das colunas
    df['Overlaped_Element_ID'] = df['Overlaped_Element_ID'].apply(lista_para_string)
    df.drop(columns=['start', 'end', 'contig'], inplace=True)


    # Adicionando tag
    for i, row in df.iterrows():
        if row['Overlaped_Element_ID'] == '':
            df.at[i, 'tag'] = 'unique'
        if row['Overlaped_Element_ID'] != '':
            df.at[i, 'tag'] = 'overlaped'

    df.to_csv(tax_file, sep="\t", index=False, header=True)
    print('DONE: Tag overlaped elements!', file = log)
    return(print('DONE: Tag overlaped elements!'))


class TagElements():
    def __init__(self, input_file, log):
        self.input_file = input_file
        self.log = log
    
    def run_tag_elemets(self):
        tag_elements(self.input_file, self.log)