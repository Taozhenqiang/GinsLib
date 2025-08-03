import subprocess
import platform

# 定位模式
mode = ["spp","spp/ins LC","spp/ins TC","ppp","ppp/ins LC","ppp/ins TC","ppd","ppd/ins LC","ppd/ins TC","ppk","ppk/ins LC","ppk/ins TC"]

# 判断操作系统
os_name = platform.system()

# 可执行文件路径与配置文件
path = ".\\..\\..\\build\\Bin\\GINSLIB.exe"
conf = "HG4930_Uav_open_20220125.conf"

# 根据操作系统调整命令
if os_name != "Windows":
    # 如果是Linux或MacOS，根据需要修改路径
        path = "./../build/Bin/GINSLIB"  # 更新命令路径

# Windows上的命令
commands = [
    [path,"-k",conf,"-p","0","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp
    [path,"-k",conf,"-p","0","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp/ins LC
    [path,"-k",conf,"-p","0","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp/ins TC
    [path,"-k",conf,"-p","8","-gins","0","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp
    [path,"-k",conf,"-p","8","-gins","1","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp/ins LC
    [path,"-k",conf,"-p","8","-gins","2","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp/ins TC
    [path,"-k",conf,"-p","1","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppd
    [path,"-k",conf,"-p","1","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppp/ins LC
    [path,"-k",conf,"-p","1","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppp/ins TC
    [path,"-k",conf,"-p","2","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppk
    [path,"-k",conf,"-p","2","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppp/ins LC
    [path,"-k",conf,"-p","2","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppp/ins TC

    # [path,"-k",conf,"-p","2","-gins","2","-ion","1","-tro","1","-eph","0","-flt","1","-amb","2"], #ppp/ins TC IGG3
] 

# 循环执行多个命令
for command, m in zip(commands, mode):
    # 执行命令
    result = subprocess.run(command)

    # 检查命令执行结果
    if result.returncode == 0:
        print(f"{m} is ok!")
    else:
        print(f"error: {result.returncode}")
