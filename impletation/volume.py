import pandas as pd
import numpy as np
import sys

# 设置编码以支持中文输出
sys.stdout.reconfigure(encoding='utf-8')


def calculate_co2_saturation(porosity):
    """
    Calculate CO2 saturation based on porosity using the formula:
    S_CO2 = -0.3136 * ln(φ) - 0.1334

    Args:
        porosity: Reservoir porosity (as a fraction, e.g., 0.15 for 15%)

    Returns:
        CO2 saturation (as a fraction)
    """
    return -0.3136 * np.log(porosity) - 0.1334


def estimate_co2_density(depth, temperature=None, pressure=None):
    """
    Estimate CO2 density based on depth, temperature, and pressure.
    If temperature and pressure are not provided, they are estimated based on depth.

    Args:
        depth: Depth in meters
        temperature: Temperature in °C (optional)
        pressure: Pressure in MPa (optional)

    Returns:
        CO2 density in kg/m³
    """
    # If temperature is not provided, estimate it based on depth
    # Assuming geothermal gradient of 25°C/km and surface temperature of 15°C
    if temperature is None:
        temperature = 15 + 0.025 * depth

    # If pressure is not provided, estimate it based on depth
    # Assuming hydrostatic pressure gradient of 10 MPa/km
    if pressure is None:
        pressure = 0.1 * depth

    # Simplified density model based on temperature and pressure
    # This is a simplified model and should be replaced with a more accurate one
    # for real-world applications
    if temperature < 31.1 and pressure > 7.38:  # Critical point of CO2
        # Supercritical CO2
        density = 600 + 70 * (pressure - 7.38) - 10 * (temperature - 31.1)
    else:
        # Gaseous CO2 (simplified)
        density = 1.98 * pressure * 1000 / (0.08314 * (temperature + 273.15))

    return max(50, min(density, 900))  # Constrain to reasonable values


def estimate_co2_solubility(temperature, pressure, salinity=0):
    """
    Estimate CO2 solubility in water based on temperature, pressure, and salinity.

    Args:
        temperature: Temperature in °C
        pressure: Pressure in MPa
        salinity: Salinity in weight percent (default: 0)

    Returns:
        CO2 solubility in mol/kg
    """
    # Simplified solubility model based on temperature and pressure
    # This is a simplified model and should be replaced with a more accurate one
    # for real-world applications

    # Base solubility at standard conditions
    base_solubility = 0.03  # mol/kg at 25°C, 0.1 MPa

    # Pressure effect (increases with pressure)
    pressure_factor = 0.01 * pressure

    # Temperature effect (decreases with temperature)
    temperature_factor = max(0, 1 - 0.01 * (temperature - 25))

    # Salinity effect (decreases with salinity)
    salinity_factor = max(0.5, 1 - 0.05 * salinity)

    solubility = base_solubility * pressure_factor * temperature_factor * salinity_factor

    return max(0.01, min(solubility, 1.0))  # Constrain to reasonable values


def calculate_bound_storage(area, thickness, porosity, co2_density=None, co2_saturation=None, depth=2000):
    """
    Calculate bound gas theoretical storage capacity.

    M_bound = A × H × φ × S_CO2 × ρ_CO2 / 10^3

    Args:
        area: Reservoir distribution area (km²)
        thickness: Average effective reservoir thickness (m)
        porosity: Average reservoir porosity (as a fraction)
        co2_density: CO2 density under reservoir conditions (kg/m³)
        co2_saturation: CO2 saturation (as a fraction)
        depth: Average reservoir depth (m)

    Returns:
        Bound gas theoretical storage capacity (Mt)
    """
    # Calculate CO2 saturation if not provided
    if co2_saturation is None:
        co2_saturation = calculate_co2_saturation(porosity)

    # Estimate CO2 density if not provided
    if co2_density is None:
        co2_density = estimate_co2_density(depth)

    # Calculate bound gas storage capacity (Mt)
    m_bound = area * thickness * porosity * co2_saturation * co2_density / 1e3

    return m_bound


def calculate_dissolved_storage(area, thickness, porosity, co2_saturation=None,
                              water_density=1000, co2_solubility=None,
                              depth=2000, temperature=None, pressure=None, salinity=3):
    """
    Calculate dissolved gas theoretical storage capacity.

    M_dissolved = A × H × φ × ρ_w × R_CO2 × M_CO2 × (1 - S_CO2) / 10^3

    Args:
        area: Reservoir distribution area (km²)
        thickness: Average effective reservoir thickness (m)
        porosity: Average reservoir porosity (as a fraction)
        co2_saturation: CO2 saturation (as a fraction)
        water_density: Formation water density (kg/m³)
        co2_solubility: CO2 solubility in formation water (mol/kg)
        depth: Average reservoir depth (m)
        temperature: Reservoir temperature (°C)
        pressure: Reservoir pressure (MPa)
        salinity: Formation water salinity (weight percent)

    Returns:
        Dissolved gas theoretical storage capacity (Mt)
    """
    # Calculate CO2 saturation if not provided
    if co2_saturation is None:
        co2_saturation = calculate_co2_saturation(porosity)

    # Estimate temperature and pressure if not provided
    if temperature is None:
        temperature = 15 + 0.025 * depth  # Geothermal gradient

    if pressure is None:
        pressure = 0.1 * depth  # Hydrostatic pressure

    # Estimate CO2 solubility if not provided
    if co2_solubility is None:
        co2_solubility = estimate_co2_solubility(temperature, pressure, salinity)

    # Calculate dissolved gas storage capacity (Mt)
    m_dissolved = area * thickness * porosity * water_density * co2_solubility * CO2_MOLAR_MASS * (1 - co2_saturation) / 1e3

    return m_dissolved


def calculate_mineral_storage(area, thickness, rock_density=2650, reactive_mineral_fraction=0.1,
                            reaction_efficiency=0.2, co2_fixation_factor=0.58):
    """
    Calculate mineral fixation theoretical storage capacity.

    M_mineral = A × H × ρ_rock × C_react × η × f_CO2

    Args:
        area: Reservoir distribution area (km²)
        thickness: Average effective reservoir thickness (m)
        rock_density: Rock density (kg/m³)
        reactive_mineral_fraction: Mass fraction of reactive minerals (%)
        reaction_efficiency: Mineral reaction efficiency (%)
        co2_fixation_factor: CO2 fixation amount per unit mineral (kg/kg)

    Returns:
        Mineral fixation theoretical storage capacity (Mt)
    """
    # Calculate mineral fixation storage capacity (Mt)
    m_mineral = area * thickness * rock_density * reactive_mineral_fraction * reaction_efficiency * co2_fixation_factor / 1e3

    return m_mineral


def calculate_effective_storage(m_bound, m_dissolved, m_mineral, efficiency_factor=0.015):
    """
    Calculate effective storage capacity.

    M_effective = (M_bound + M_dissolved + M_mineral) × C_e

    Args:
        m_bound: Bound gas theoretical storage capacity (Mt)
        m_dissolved: Dissolved gas theoretical storage capacity (Mt)
        m_mineral: Mineral fixation theoretical storage capacity (Mt)
        efficiency_factor: Efficiency coefficient (0.01-0.02 for basin level)

    Returns:
        Effective storage capacity (Mt)
    """
    # Calculate effective storage capacity (Mt)
    m_effective = (m_bound + m_dissolved + m_mineral) * efficiency_factor

    return m_effective


def load_config():
    """
    从配置文件中加载参数。

    Returns:
        包含配置参数的字典
    """
    try:
        config_df = pd.read_csv('data/config.csv')
        config = {}

        # 参数名称映射，将中文参数名映射到英文参数名
        param_mapping = {
            '有效系数': 'efficiency_factor',
            '二氧化碳摩尔质量': 'co2_molar_mass'
        }

        for _, row in config_df.iterrows():
            param_name = row['参数']
            if param_name in param_mapping:
                config[param_mapping[param_name]] = row['数值']

        return config
    except Exception as e:
        print(f"加载配置文件时出错: {e}")
        # 如果无法加载配置文件，使用默认值
        return {
            'efficiency_factor': 0.015,
            'co2_molar_mass': 0.044
        }


# 定义全局常量
CO2_MOLAR_MASS = 0.044  # kg/mol, 默认值
EFFICIENCY_FACTOR = 0.015  # 默认值

# 加载配置并初始化常量
def _initialize_constants():
    """从配置文件加载并更新全局常量"""
    config = load_config()
    global CO2_MOLAR_MASS, EFFICIENCY_FACTOR
    CO2_MOLAR_MASS = float(config.get('co2_molar_mass', 0.044))  # kg/mol
    EFFICIENCY_FACTOR = float(config.get('efficiency_factor', 0.015))

# 初始化常量
_initialize_constants()


def load_basin_data(file_path='data/china_saline_aquifers.csv'):
    """
    从 CSV 文件中加载盆地数据。

    Args:
        file_path: 包含盆地数据的 CSV 文件路径

    Returns:
        包含盆地数据的 DataFrame
    """
    try:
        df = pd.read_csv(file_path)

        # 列名映射，将中文列名映射到英文列名
        column_mapping = {
            '盆地名称': 'basin_name',
            '面积（平方公里）': 'area_km2',
            '厚度（米）': 'thickness_m',
            '孔隙度': 'porosity',
            '深度（米）': 'depth_m',
            '盐度（重量百分比）': 'salinity_wt_percent',
            '可反应矿物质量分数': 'reactive_mineral_fraction',
            '矿物反应效率': 'reaction_efficiency',
            '岩石密度（千克每立方米）': 'rock_density_kg_m3',
            '地层水密度（千克每立方米）': 'water_density_kg_m3',
            '二氧化碳固定系数': 'co2_fixation_factor'
        }

        # 重命名列
        df = df.rename(columns=column_mapping)

        return df
    except Exception as e:
        print(f"加载盆地数据时出错: {e}")
        return None


def estimate_china_saline_aquifer_storage(data_file_path='data/china_saline_aquifers.csv'):
    """
    估算中国深层咸水层的二氧化碳总存储容量。

    Args:
        data_file_path: 包含盆地数据的 CSV 文件路径

    Returns:
        包含不同盆地存储容量及总量的字典
    """
    # 从 CSV 文件加载盆地数据
    basin_df = load_basin_data(data_file_path)

    if basin_df is None:
        print("无法加载盆地数据。使用默认值。")
        return {}

    # 计算每个盆地的存储容量
    results = {}
    total_bound = 0
    total_dissolved = 0
    total_mineral = 0
    total_effective = 0

    for _, row in basin_df.iterrows():
        # Extract parameters from the dataframe row
        basin_name = row['basin_name']
        area = row['area_km2']
        thickness = row['thickness_m']
        porosity = row['porosity']
        depth = row['depth_m']
        salinity = row['salinity_wt_percent']
        reactive_mineral_fraction = row['reactive_mineral_fraction']
        reaction_efficiency = row['reaction_efficiency']
        rock_density = row['rock_density_kg_m3']
        water_density = row['water_density_kg_m3']
        co2_fixation_factor = row['co2_fixation_factor']

        # Estimate temperature and pressure
        temperature = 15 + 0.025 * depth
        pressure = 0.1 * depth

        # Calculate CO2 saturation
        co2_saturation = calculate_co2_saturation(porosity)

        # Estimate CO2 density
        co2_density = estimate_co2_density(depth, temperature, pressure)

        # Estimate CO2 solubility
        co2_solubility = estimate_co2_solubility(temperature, pressure, salinity)

        # Calculate storage capacities
        m_bound = calculate_bound_storage(
            area, thickness, porosity, co2_density, co2_saturation, depth
        )

        m_dissolved = calculate_dissolved_storage(
            area, thickness, porosity, co2_saturation, water_density, co2_solubility,
            depth, temperature, pressure, salinity
        )

        m_mineral = calculate_mineral_storage(
            area, thickness, rock_density, reactive_mineral_fraction, reaction_efficiency, co2_fixation_factor
        )

        m_effective = calculate_effective_storage(m_bound, m_dissolved, m_mineral, EFFICIENCY_FACTOR)

        # Store results
        results[basin_name] = {
            "bound_storage": m_bound,
            "dissolved_storage": m_dissolved,
            "mineral_storage": m_mineral,
            "effective_storage": m_effective,
            "parameters": {
                "area": area,
                "thickness": thickness,
                "porosity": porosity,
                "depth": depth,
                "co2_saturation": co2_saturation,
                "co2_density": co2_density,
                "co2_solubility": co2_solubility
            }
        }

        # Add to totals
        total_bound += m_bound
        total_dissolved += m_dissolved
        total_mineral += m_mineral
        total_effective += m_effective

    # Add totals to results
    results["Total"] = {
        "bound_storage": total_bound,
        "dissolved_storage": total_dissolved,
        "mineral_storage": total_mineral,
        "effective_storage": total_effective
    }

    return results


def print_storage_results(results):
    """
    以格式化表格的形式打印存储容量结果。

    Args:
        results: 包含存储容量结果的字典
    """
    print("\n中国深层咸水层二氧化碳存储容量估算\n")
    print("{:<25} {:<15} {:<15} {:<15} {:<15}".format(
        "盆地", "束缚气（Mt）", "溶解气（Mt）", "矿物固化（Mt）", "有效存储（Mt）"
    ))
    print("-" * 85)

    for basin_name, data in results.items():
        if basin_name != "Total":
            print("{:<25} {:<15.2f} {:<15.2f} {:<15.2f} {:<15.2f}".format(
                basin_name,
                data["bound_storage"],
                data["dissolved_storage"],
                data["mineral_storage"],
                data["effective_storage"]
            ))

    print("-" * 85)
    print("{:<25} {:<15.2f} {:<15.2f} {:<15.2f} {:<15.2f}".format(
        "总计",
        results["Total"]["bound_storage"],
        results["Total"]["dissolved_storage"],
        results["Total"]["mineral_storage"],
        results["Total"]["effective_storage"]
    ))


def export_results_to_csv(results, output_file='data/co2_storage_results.csv'):
    """
    将存储容量结果导出到 CSV 文件。

    Args:
        results: 包含存储容量结果的字典
        output_file: 输出 CSV 文件的路径
    """
    # 从结果创建 DataFrame
    data = []
    for basin_name, basin_data in results.items():
        if basin_name != "Total":
            data.append({
                '盆地': basin_name,
                '束缚气存储量（Mt）': basin_data['bound_storage'],
                '溶解气存储量（Mt）': basin_data['dissolved_storage'],
                '矿物固化存储量（Mt）': basin_data['mineral_storage'],
                '有效存储量（Mt）': basin_data['effective_storage']
            })

    # 添加总计行
    data.append({
        '盆地': '总计',
        '束缚气存储量（Mt）': results['Total']['bound_storage'],
        '溶解气存储量（Mt）': results['Total']['dissolved_storage'],
        '矿物固化存储量（Mt）': results['Total']['mineral_storage'],
        '有效存储量（Mt）': results['Total']['effective_storage']
    })

    # 创建并保存 DataFrame
    df = pd.DataFrame(data)
    try:
        df.to_csv(output_file, index=False)
        print(f"\n结果已导出到 {output_file}")
        return True
    except Exception as e:
        print(f"\n导出结果到 CSV 文件时出错: {e}")
        return False


def main(data_file_path=None, export_csv=True, output_file='data/co2_storage_results.csv'):
    """
    估算和显示中国深层咸水层二氧化碳存储容量的主函数。

    Args:
        data_file_path: 包含盆地数据的 CSV 文件路径（可选）
        export_csv: 是否将结果导出到 CSV 文件（默认：True）
        output_file: 输出 CSV 文件的路径（默认：'data/co2_storage_results.csv'）

    Returns:
        包含存储容量结果的字典
    """
    print("正在估算中国深层咸水层二氧化碳存储容量...")

    # 估算存储容量
    if data_file_path:
        print(f"使用数据来源: {data_file_path}")
        results = estimate_china_saline_aquifer_storage(data_file_path)
    else:
        print(f"使用默认数据来源: data/china_saline_aquifers.csv")
        results = estimate_china_saline_aquifer_storage()

    # 打印结果
    print_storage_results(results)

    # 如果需要，将结果导出到 CSV 文件
    if export_csv and results:
        export_results_to_csv(results, output_file)

    print("\n注意: 这些是基于简化模型和示例数据的估计值。")
    print("要进行准确评估，需要详细的地质数据和更复杂的模型。")

    # 返回结果以便进一步分析
    return results


if __name__ == "__main__":
    main()
