### 对getrusage, time, /proc/pid/status三种统计RAM方法的分析

#### execvp + getrusage
    [bwamem-main] Real time: 180.227 sec; CPU: 3442.331 sec
        Max_RAM: 7617864 kbytes
        Max_VIRT:  8165064 kbytes
        MAX_rss: 7617864 kbytes
        CPU Seconds: 3442.351
    [toy-main] Real time: 180.253 sec; CPU: 83.077 sec
    
#### /usr/bin/time -v  
    [main] Real time: 219.194 sec; CPU: 4155.078 sec
        User time (seconds): 4128.21
        System time (seconds): 26.88
        Elapsed (wall clock) time (h:mm:ss or m:ss): 3:39.21
        Maximum resident set size (kbytes): 30474992

```
time和toy pstat的CPU时间和时钟时间都与程序消耗几乎一致.
```

Table: RAM usage of 5 mappers on ERP001775_s4_1

|Mapper|status|getrusage|time|
|------|------|---------|----|
|bowtie2|5,736|3,569,788|13,972,192|
|bwamem|7,957,828|7,957,828|31,721,936|
|gem|14,056,732|14,056,736|56,247,504|
|kart|11,546,256|11,546,256|46,185,296|
|ZM-S|57,089,712|57,089,716|228,235,904|
|ZM-M|28,627,400|28,627,400|114,873,264|
|ZM-P|14,708,336|14,708,336|59,356,128|

    系统调用getrusage与手动访问(100us每次)/proc/pid/status文件中的"VmRSS"最大值几乎一致.
    奇怪的是bowtie2的手动访问结果明显是错误的, 具体原因未知.
    getrusage与top命令的查看情况几乎吻合, 可以被信任.
    time给出的结果始终偏大, 将其除以getrusage得到: 3.914, 3.986, 4.001, 4.000, 3.998, 4.013, 4.036
    time输出的MRSS大约是其4倍, 据网络资料所言, 这很可能是程序bug.
    


### Term

#### RSS
    RSS: Resident set size is the amount of memory the residents(presents) in (real) RAM.
    MRSS: Maximum RSS is peak RAM memory usage.
    有趣的是, linux上所有进程的RSS加并不等于这些程序实际的资源使用,
    因为RSS多次记录了shared libraries, 所以大于等于实际资源消耗.
    
#### VIRT
    VIRT >= real RAM + swapped memory.