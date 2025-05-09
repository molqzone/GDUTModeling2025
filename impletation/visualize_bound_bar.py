import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib

# 设置中文字体，自动适配常见操作系统
zh_fonts = ["SimHei", "Microsoft YaHei", "STSong", "Arial Unicode MS"]
for font in zh_fonts:
    if font in matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf') or font in [f.name for f in matplotlib.font_manager.fontManager.ttflist]:
        plt.rcParams['font.sans-serif'] = [font]
        break
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(script_dir, '../data/china_saline_aquifers_with_bound.csv')

# 读取数据，跳过格式错误的行
try:
    df = pd.read_csv(data_path, encoding='utf-8-sig', on_bad_lines='skip')
except TypeError:
    # 兼容旧版pandas
    df = pd.read_csv(data_path, encoding='utf-8-sig', error_bad_lines=False)

# 排除“总计”行
df = df[df['地区'] != '总计']

# 只保留有储存量的地区
df = df.dropna(subset=['束缚气储存量_吨'])

# 按储存量降序排序
df = df.sort_values('束缚气储存量_吨', ascending=False)

plt.figure(figsize=(12, 6))
bars = plt.bar(df['地区'], df['束缚气储存量_吨'] / 1e9)  # 单位转为亿吨
plt.ylabel('束缚气储存量 (亿吨)')
plt.xlabel('地区')
plt.title('中国各地区咸水层束缚气储存量柱状图')
plt.xticks(rotation=45, ha='right')

# 在柱子顶部显示数据文本
for bar, value in zip(bars, df['束缚气储存量_吨'] / 1e9):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{value:.2f}',
             ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.show()
