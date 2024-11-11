# COG_category sütunundaki 'AB' değerini ayırarak listeye çevirme
df['COG_category'] = df['COG_category'].apply(lambda x: list(x) if len(x) > 1 else [x])
# explode fonksiyonunu kullanarak her kategoriye ayrı satır olarak ayırma
df_exploded = df.explode('COG_category')
# Sadece '#query' ve 'COG_category' sütunlarını seçme
df_COG = df_exploded[['#query', 'COG_category']]
df_COG=df_COG[df_COG['COG_category'] != "-"]
# Sonucu görüntüleme
print(df_COG)


#df.explode('COG_category')
#df['COG_category'] = df['COG_category'].str.split(',')
#df = df.explode('COG_category')
#df = df[df['COG_category'] != "-"]
