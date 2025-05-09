import pandas as pd
import numpy as np
import os

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(script_dir, '../data/china_saline_aquifers.csv')

# 读取CSV数据（新表头）
df = pd.read_csv(data_path, encoding='utf-8')

# 打印字段名调试
print('字段名:', list(df.columns))

# 束缚气储存量计算公式
# M_bound = A * H * phi * S_CO2 * rho_CO2 / 1e3
# rho_CO2: 取自表格（CO2密度/kg/m³），如缺失则默认725

def safe_float(val, default=np.nan):
    try:
        if pd.isna(val):
            return default
        if isinstance(val, str):
            val = val.strip().replace('\u3000', '').replace(' ', '')  # 去除全角/半角空格
            if val == '' or val.lower() == 'nan':
                return default
            if '-' in val:  # 处理区间
                val = val.split('-')[0]
        return float(val)
    except Exception:
        return default

def calc_M_bound(row):
    # 优先使用深部含水层储量
    V = safe_float(row.get('深部含水层储量/10^8m³'))  # 10^8 m³
    phi = safe_float(row.get('孔隙度/%'), 15)  # %
    rho_CO2 = safe_float(row.get('CO2密度/kg/m³'), 725)  # kg/m³
    S_CO2 = safe_float(row.get('CO2饱和度/%'))
    if not np.isnan(V):
        phi_frac = phi / 100
        S_CO2_frac = S_CO2
        if np.isnan(S_CO2_frac):
            S_CO2_frac = 0.4
        M_bound = V * 1e8 * phi_frac * S_CO2_frac * rho_CO2 / 1e3
        return M_bound
    # 否则用面积和厚度
    A = safe_float(row.get('地下含盐水面积/km²'))  # km²
    H = safe_float(row.get('厚度/m'))  # m
    if np.isnan(A) or np.isnan(H):
        return np.nan
    phi_frac = phi / 100
    S_CO2_frac = S_CO2
    if np.isnan(S_CO2_frac):
        S_CO2_frac = 0.4
    M_bound = A * 1e6 * H * phi_frac * S_CO2_frac * rho_CO2 / 1e3
    return M_bound

df['束缚气储存量_吨'] = df.apply(calc_M_bound, axis=1)

# 输出结果
print(df[['地区', '束缚气储存量_吨']])
print('总计:', df['束缚气储存量_吨'].sum())

# 保存结果并添加总计行
output_path = os.path.join(script_dir, '../data/china_saline_aquifers_with_bound.csv')
df.to_csv(output_path, index=False, encoding='utf-8-sig')

# 追加总计行到csv，字段数与表头一致
with open(output_path, 'a', encoding='utf-8-sig') as f:
    total_cols = len(df.columns)
    f.write('总计' + ',' * (total_cols - 2) + f',{df["束缚气储存量_吨"].sum()}' + '\n')
