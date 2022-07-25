import pandas as pd
import re
import os

def mergeCSV(data_path, qtd_skip_rows):
    dfO = pd.DataFrame()
    
    fList = os.listdir(r'' + data_path + '')
    for f in fList:
        path = data_path + "/" + f
        
        df = pd.read_csv((path), skiprows=qtd_skip_rows)
        country_substring = re.search(r'(.+?)-', f)
        country = country_substring.group(1)
        df['Country'] = country
        if dfO.empty:
            dfO = df.copy()
        else:
            frames = [dfO, df]
            dfO =  pd.concat(frames)
    dfO.to_csv(data_path + '_all.csv', encoding='utf-8')
    
mergeCSV('gdp', 16)
mergeCSV('healthcare', 16)
mergeCSV('population', 15)