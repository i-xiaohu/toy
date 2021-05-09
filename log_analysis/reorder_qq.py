from log_analysis.utils import *


def fmt2sec(f):
    x = f.split(":")
    ret = int(x[0]) * 3600
    ret += int(x[1]) * 60
    ret += int(x[2])
    return ret


def sec2fmt(s):
    s = math.ceil(float(s))
    x = [0, 0, 0]
    x[0] = s // 3600
    s -= x[0] * 3600
    x[1] = s // 60
    s -= x[1] * 60
    x[2] = s
    return "{:d}:{:02d}:{:02d}".format(x[0], x[1], x[2])


# Format RAM from 'pstat'
def fmt_ram(var):
    if var == 0: return "None"
    var = float(var)
    if var > 1024 * 1024:
        return "{:d}G".format(math.ceil(var/1024/1024))
    elif var > 1024:
        return "{:d}M".format(math.ceil(var/1024))
    else:
        return "<1M"


if __name__ == "__main__":
    datasets = ['human_E1_L125_25X_1', 'human_E1_L150_25X_1']
    os.chdir(r'D:\ClionProject\zsmem-experiments\log\toy_sort')
    table = []
    table.append("|Dataset|Comp|Reads|Input|Trie|Search|Write|Matched|Same|RAM|Total|")
    table.append("|-------|----|-----|-----|----|------|-----|-------|----|---|-----|")
    for ds in datasets:
        for fn in os.listdir("."):
            if fn.find(ds) == -1: continue
            f = open(fn, "r")
            reads_n, reading_time = 0, 0
            match_reads, chunk_time, search_time = 0, 0, 0
            same_reads, same_percent = 0, 0
            ram, total = 0, 0
            while True:
                s = f.readline()
                if s == '': break
                sp = s.strip().split()
                if s.find("Input") != -1:
                    reads_n = sp[1]
                    reading_time = fmt2sec(sp[-2])
                if s.find("Chunk") != -1:
                    chunk_time += fmt2sec(sp[-2])
                if s.find("read matched") != -1:
                    match_reads += int(sp[0].replace(",", ""))
                    search_time += fmt2sec(sp[-2])
                if s.find("same reads") != -1:
                    same_reads = sp[1]
                    same_percent = sp[5]
                if s.find("MAX_rss") != -1:
                    ram = fmt_ram(sp[-2])
                if s.find("Elapsed_seconds") != -1:
                    total = math.ceil(float(sp[-1]))
            write = total - reading_time - chunk_time - search_time
            table.append("|{ds}|{comp}|{reads}|{input}|{trie}|{search}|{write}|{match}|{same}|{ram}|{total}|".format(
                ds=ds,
                comp=fn.split('.')[0],
                reads=reads_n,
                input=sec2fmt(reading_time),
                trie=sec2fmt(chunk_time),
                search=sec2fmt(search_time),
                write=sec2fmt(write),
                match="{:,}".format(match_reads),
                same=same_reads,
                ram=ram,
                total=sec2fmt(total)
            ))
            table.append("| | | |{input:.2f}%|{trie:.2f}%|{search:.2f}%|{write:.2f}%|{match:.2f}%|{same}| | |".format(
                input=100.0*reading_time/total,
                trie=100.0*chunk_time/total,
                search=100.0*search_time/total,
                write=100.0*write/total,
                match=100.0*match_reads/int(reads_n.replace(",", "")),
                same=same_percent
            ))
            put_table(r'D:\ClionProject\toy\dev_notes\sort_reads.md', "Table: Profiling", table)
            # print(reads_n, reading_time)
            # print("{:,} reads matched".format(match_reads))
            # print(sec2fmt(chunk_time), sec2fmt(search_time))
            # print(same_reads, same_percent)
            # print(ram, time)
