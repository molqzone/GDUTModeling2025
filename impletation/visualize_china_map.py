import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 支持中文
matplotlib.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 读取数据
script_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(script_dir, '../data/china_saline_aquifers_with_bound.csv')
# 跳过解析错误的行（如总计行）
df = pd.read_csv(data_path, encoding='utf-8', on_bad_lines='skip')
df = df[~df['地区'].astype(str).str.contains('总计', na=False)]

# 只保留有储存量的地区
df = df.dropna(subset=['束缚气储存量_吨'])

# 读取地区-省级行政区映射表
map_config_path = os.path.join(script_dir, '../data/map_config.csv')
map_df = pd.read_csv(map_config_path, encoding='utf-8')
df = df.merge(map_df, on='地区', how='left')

# 省级行政区地理数据（可用 geopandas 自带的或下载中国省级 shapefile）
# 这里假设有一个 gadm41_CHN_1.shp 文件在 data 目录下
shp_path = os.path.join(script_dir, '../data/gadm41_CHN_shp/gadm41_CHN_1.shp')
gdf = gpd.read_file(shp_path)

# 合并数据，使用标准省份名匹配
merged = gdf.merge(df, left_on='NAME_1', right_on='省级行政区', how='left')

# 输出合并后的主要字段和匹配成功的省份数
print(merged[['NAME_1', '地区', '束缚气储存量_吨']])
print('匹配成功省份数：', merged['束缚气储存量_吨'].notna().sum())

# 绘制地图
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
merged.plot(column='束缚气储存量_吨', ax=ax, legend=True, cmap=cm.OrRd, missing_kwds={
    "color": "lightgrey",
    "label": "No data"
})
ax.set_title('中国各地区咸水层CO$_2$束缚气储存量分布', fontsize=16)
ax.axis('off')
plt.tight_layout()
# 在有数据的省份中心添加省份名和储存量文本
for idx, row in merged.iterrows():
    if pd.notna(row['束缚气储存量_吨']):
        x, y = row['geometry'].centroid.x, row['geometry'].centroid.y
        label = f"{row['NAME_1']}\n{row['束缚气储存量_吨']:.2e}"
        ax.text(x, y, label, fontsize=8, ha='center', va='center', color='black')
plt.show()
