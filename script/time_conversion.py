import math
from datetime import datetime, timedelta

# GPS时间参考点 (1980年1月6日 00:00:00)
gpst0 = [1980, 1, 6, 0, 0, 0]

# 闰秒表 (年份, 月份, 日, 时, 分, 秒, UTC-GPST差值)
leaps = [
    [2017, 1, 1, 0, 0, 0, -18],
    [2015, 7, 1, 0, 0, 0, -17],
    [2012, 7, 1, 0, 0, 0, -16],
    [2009, 1, 1, 0, 0, 0, -15],
    [2006, 1, 1, 0, 0, 0, -14],
    [1999, 1, 1, 0, 0, 0, -13],
    [1997, 7, 1, 0, 0, 0, -12],
    [1996, 1, 1, 0, 0, 0, -11],
    [1994, 7, 1, 0, 0, 0, -10],
    [1993, 7, 1, 0, 0, 0, -9],
    [1992, 7, 1, 0, 0, 0, -8],
    [1991, 1, 1, 0, 0, 0, -7],
    [1990, 1, 1, 0, 0, 0, -6],
    [1988, 1, 1, 0, 0, 0, -5],
    [1985, 7, 1, 0, 0, 0, -4],
    [1983, 7, 1, 0, 0, 0, -3],
    [1982, 7, 1, 0, 0, 0, -2],
    [1981, 7, 1, 0, 0, 0, -1],
    [0, 0, 0, 0, 0, 0, 0]  # 结束标记
]

class GTimeT:
    """Python版本的gtime_t结构体"""
    def __init__(self, time=None, sec=0.0):
        self.time = time  # 时间戳 (秒)
        self.sec = sec    # 秒的小数部分
        
    def __repr__(self):
        return f"GTimeT(time={self.time}, sec={self.sec})"

def epoch2time(ep):
    """
    将日历日期时间转换为GTimeT结构体
    args: ep - [年, 月, 日, 时, 分, 秒]
    return: GTimeT对象
    """
    # 强制转换为整数类型，与C语言版本保持一致
    year = int(ep[0])
    mon = int(ep[1])
    day = int(ep[2])
    hour = int(ep[3])
    minute = int(ep[4])
    sec_val = ep[5]  # 秒可以是浮点数
    
    # 检查年份范围 (1970-2099)
    if year < 1970 or year > 2099 or mon < 1 or mon > 12:
        return GTimeT(0, 0.0)
    
    # 计算从1970年开始的天数
    doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    days = (year - 1970) * 365 + (year - 1969) // 4 + doy[mon - 1] + day - 2
    # 闰年调整 (1901-2099年间能被4整除的年份为闰年)
    if year % 4 == 0 and mon >= 3:
        days += 1
    
    sec_int = int(math.floor(sec_val))  # 秒的整数部分
    total_seconds = days * 86400 + hour * 3600 + minute * 60 + sec_int
    sec_frac = sec_val - sec_int
    
    return GTimeT(total_seconds, sec_frac)

def timeadd(t, sec):
    """
    时间加法运算
    args: t - GTimeT对象, sec - 要加的秒数
    return: GTimeT对象 (t + sec)
    """
    result = GTimeT(t.time, t.sec + sec)
    tt = math.floor(result.sec)
    result.time += int(tt)
    result.sec -= tt
    return result

def timediff(t1, t2):
    """
    时间差运算
    args: t1, t2 - GTimeT对象
    return: 时间差 (t1 - t2) 秒
    """
    return (t1.time - t2.time) + (t1.sec - t2.sec)

def time2gpst(t, week=None):
    """
    将GTimeT结构体转换为GPS周和周内秒
    args: t - GTimeT对象, week - 周数输出变量 (可为None)
    return: GPS周内秒
    """
    t0 = epoch2time(gpst0)
    sec = t.time - t0.time
    w = int(sec // (86400 * 7))
    
    if week is not None:
        week[0] = w
    
    return (sec - w * 86400 * 7) + t.sec

def utc2gpst(t):
    """
    将UTC时间转换为GPS时间，考虑闰秒
    args: t - UTC时间的GTimeT对象
    return: GPS时间的GTimeT对象
    """
    for leap in leaps:
        if leap[0] == 0:  # 结束标记
            break
        
        leap_time = epoch2time(leap[:6])
        if timediff(t, leap_time) >= 0.0:
            return timeadd(t, -leap[6])
    
    return t