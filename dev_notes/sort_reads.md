排序目标: 784,363,613条序列; 78G碱基.  

    程序总共运行23:52:14  
    输入目标花费00:16:05  
    建树和查询共20:21:47  
    目前以1G作为构建trie的字符上限, 总共运行79轮, 每轮trie大约占68G, reads占约2G, 平均时间为00:15:33  
    到最后几轮时, 建树和查询的时间大约为5分钟.
    写入(gzip pipied)的消耗为3:14:22
    程序内存消耗峰值为328G, 属于可以接受的范畴.
    
可见建树的时间消耗并不大, 主要是查询消耗, 需要并行支持来对查询进行加速.  

Found 2,405,357 same reads, taking 0.31% of the reads in trie.  
相同序列属于小概率事件, 对全局解的影响不大.  