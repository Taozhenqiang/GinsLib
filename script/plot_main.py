# main.py
import os
from read_sol import read_solution
from plot_trajectory import plot_trajectory
from plot_position   import plot_position
from plot_velocity   import plot_velocity
from plot_attitude   import plot_attitude
from plot_imubias    import plot_imubias

def main():
    # 提供文件路径和文件名
    pathname = 'G:/Navigation_Learn/GNSS/RTKLIB-b34k/data/STIM300_Vehicle_complex_20230102/result/'
    filename = 'STIM300_Vehicle_complex_20230102_PPK_TC.pos'
    navfile = os.path.join(pathname, filename)

    # 设置文件读取配置
    skip_lines = 28
    row = 10000
    # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
    col = 18
    sample =5e-2

    # 调用 read_solution 函数来读取文件数据
    data = read_solution(navfile, skip_lines, row, col, sample)

    # 检查是否成功读取数据
    if data is not None:
        print("文件数据读取成功")
        print(data)
    else:
        print("文件数据读取失败")
        return

    # 调用 plot_trajectory 函数来绘制轨迹
    # plot_trajectory(data)

    # plot_position(data)

    # plot_velocity(data)

    # plot_attitude(data)

    plot_imubias(data)

if __name__ == "__main__":
    main()
