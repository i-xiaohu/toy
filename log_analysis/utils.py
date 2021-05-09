import math
import os


# Get accuracy percentage.
def get_accuracy(fn):
    if os.path.exists(fn) is False: return 0
    f = open(fn, "r")
    if f is None:
        return 0
    unmap = int(f.readline().strip())
    correct, mapped = 0, 0
    while True:
        line = f.readline().strip()
        if line == '':
            break
        mapped += int(line.split()[1])
        correct += int(line.split()[2])
    reads_n = mapped + unmap
    f.close()
    return 100.0 * correct / reads_n


# Get correctly aligned reads.
def get_correct(fn):
    if os.path.exists(fn) is False: return 0
    f = open(fn, "r")
    if f is None:
        return 0
    f.readline()
    correct = 0
    while True:
        line = f.readline().strip()
        if line == '':
            break
        correct += int(line.split()[2])
    f.close()
    return correct


# Get correctly aligned reads.
def get_mapped(fn):
    if os.path.exists(fn) is False: return 0
    f = open(fn, "r")
    if f is None:
        return 0
    f.readline()
    mapped = 0
    while True:
        line = f.readline().strip()
        if line == '':
            break
        mapped += int(line.split()[1])
    f.close()
    return mapped


# Return #reads of the dataset
def get_reads(ds):
    f = open("../../results/dataset.md", "r", encoding='utf-8')
    while True:
        s = f.readline()
        if s == "": break
        if s.find("Table: NGS datasets") != -1:
            f.readline()
            break
    header = [_.strip() for _ in f.readline().strip().split("|")]
    f.readline()
    while True:
        s = f.readline()
        if s == "" or s == "\n": break
        row = [_.strip() for _ in s.strip().split("|")]
        if row[header.index("ID")].find(ds) != -1:
            return int(row[header.index("#reads")].replace(',', ''))
    return -1


# Get Elapsed time from 'usr/bin/time -v'
def get_realtime(lines):
    ret = None
    for s in lines:
        if s.find("Elapsed (wall clock) time (h:mm:ss or m:ss):") != -1:
            ret = s.strip().split()[-1]
    if ret is None: return None
    return fmt_time(get_seconds(ret))


# Format seconds from 'pstat'
def fmt_time(var):
    var = math.ceil(float(var))
    hour = var // 3600
    var -= hour * 3600
    min = var // 60
    var -= min * 60
    sec = var
    return "{:d}:{:02d}:{:02d}".format(hour, min, sec)


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


# Put the table into markdown file.
def put_table(md_filename, table_name, table):
    f = open(md_filename, mode="r", encoding='utf-8')
    start = -1
    lines = []
    while True:
        s = f.readline()
        if s == '': break
        lines.append(s)
    for i in range(0, len(lines)):
        if lines[i] == "{}\n".format(table_name):
            start = i
            break
    if start == -1:
        print("'{}' not found, new table will be appended.".format(table_name))
        start = len(lines)
    end = start + 1
    for i in range(start+2, len(lines)):  # Skip a blank line.
        if lines[i].count('|') > 2:
            end = i + 1
        else:
            break
    f.close()
    if end == start + 1:
        print("'{}' is not complete, new rows will be appended.".format(table_name))

    # Output the front.
    f = open(md_filename, mode="w", encoding='utf-8')
    for i in range(0, start):
        f.write(lines[i])
    # Insert the table.
    f.write(table_name + "\n")
    f.write("\n")
    for t in table:
        f.write(t + "\n")
    # Output the tail.
    for i in range(end, len(lines)):
        f.write(lines[i])
    f.close()


# Buffer lines from the file.
def get_lines(fn):
    if os.path.exists(fn) is False: return []
    f = open(fn, "r", encoding='utf-8')
    lines = []
    while True:
        s = f.readline()
        if s == '': break
        lines.append(s)
    f.close()
    return lines


# Fetch a variable from the buffered lines.
def fetch_var(lines, var):
    for s in lines:
        sp = s.strip().split()
        if len(sp) < 2: continue
        if sp[0].find(var) != -1:
            return sp[1]
    return 0


# Get seconds from formatted time.
def get_seconds(fmt_time):
    ta = fmt_time.split(sep=':')
    ret = 0
    base = 1
    for i in range(len(ta)-1, -1, -1):
        ret += base * float(ta[i])
        base *= 60
    return ret


# Return matlab color from (R, G, B)
def matlab_color(R, G, B):
    ret = "#"
    if len(hex(R)[2:]) < 2: ret += "0"
    ret += hex(R)[2:]
    if len(hex(G)[2:]) < 2: ret += "0"
    ret += hex(G)[2:]
    if len(hex(B)[2:]) < 2: ret += "0"
    ret += hex(B)[2:]
    return ret


if __name__ == "__main__":
    print("Utils is a tools package.")
