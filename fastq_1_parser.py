def parse_fastq(filename):
    count = 0
    total_len = 0
    min_len = None
    max_len = None
    with open(filename, "r", encoding="utf-8") as f:
        while True:
            header = f.readline()

            # конец файла
            if not header:
                break

            sequence = f.readline().strip()
            plus = f.readline()
            quality = f.readline().strip()
            cur_len = len(sequence)
            count += 1
            total_len += cur_len

            #обновляем min max
            if min_len is None or cur_len < min_len:
                min_len = cur_len
            if max_len is None or cur_len > max_len:
                max_len = cur_len

    avg_len = round(total_len / count)

    return count, min_len, avg_len, max_len


if __name__ == "__main__":

    filename = "reads.fastq"
    count, min_len, avg_len, max_len = parse_fastq(filename)

    print("Общее число прочтений:", count)
    print("Минимальная длина:", min_len)
    print("Средняя длина:", avg_len)
    print("Максимальная длина:", max_len)
