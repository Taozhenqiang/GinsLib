# main.py
import os
from read_sol import read_solution
from read_ref import read_ref
from plot_trajectory import plot_trajectory
from plot_position   import plot_position
from plot_velocity   import plot_velocity
from plot_attitude   import plot_attitude
from plot_imubias    import plot_imubias
from plot_err    import plot_err

def main():
    # 提供文件路径和结果文件名
    pathname = 'G:/Navigation_Learn/GNSS/RTKLIB-b34k/data/IGM-A1_Vehicle_open_20211217-nsh/result'
    filename = 'IGM-A1_Vehicle_open_20211217-nsh_PPK_TC_PAR.pos'
    solfile = os.path.join(pathname, filename)

    # 提供文件路径和参考文件名
    pathname = 'G:/Navigation_Learn/GNSS/RTKLIB-b34k/data/IGM-A1_Vehicle_open_20211217-nsh'
    filename = 'IGMA1_open_2GNSS.txt'
    refile = os.path.join(pathname, filename)

    # 结果文件读取配置
    skip_lines = 28
    row = 10000
    # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
    col = 18
    sample = 1/200

    # 读取结果文件
    data = read_solution(solfile, skip_lines, row, col, sample)

    # 检查是否成功读取
    if data is not None:
        print("结果文件数据读取成功")
    else:
        print("结果文件数据读取失败")
        return

    # 参考文件读取配置
    skip_lines = 20
    row = 10000
    # gps week, sow, pos[x/y/z], vel[x/y/z], att[pitch/roll/heading]
    col = 11
    ins_flag = 0
    ins_interval = 1/200
    GNSS_interval = 1
    # IE-tzq
    idx = [0,1,2,3,4,8,9,10,14,15,16]
    # IE-zf
    # idx = [0,1,9,10,11,15,16,17,22,23,21]

    # 读取参考文件
    ref_data = read_ref(refile, skip_lines, row, col, ins_flag, ins_interval, GNSS_interval, idx)

    # 检查是否成功读取
    if ref_data is not None:
        print("参考文件数据读取成功")
    else:
        print("参考文件数据读取失败")
        return
    
    # 绘制轨迹
    # plot_trajectory(data)
    # 绘制位置
    # plot_position(data)
    # 绘制速度
    # plot_velocity(data)
    # 绘制姿态
    # plot_attitude(data)
    # 绘制IMU bias
    # plot_imubias(data)

    plot_err(data, ref_data, "p")

if __name__ == "__main__":
    main()
