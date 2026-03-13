def calc_gc(filename):
    gc_count = 0

    # общее число всех нуклеотидов
    total_length = 0

    with open(filename, "r", encoding="utf-8") as f:
        while True:
            header = f.readline()
            if not header:
                break
            #парсим
            sequence = f.readline().strip()
            plus = f.readline()
            quality = f.readline().strip()

            # считаем количество G и C в текущем прочтении и увелчиваем общую длину посл-ти
            gc_count += sequence.count("G") + sequence.count("C")
            total_length += len(sequence)

    if total_length == 0:
        raise ValueError("Файл пустой или не содержит последовательностей")

    gc_percent = (gc_count / total_length) * 100
    return round(gc_percent, 2)


if __name__ == "__main__":
    fname = "reads.fastq"

    gc_percent = calc_gc(fname)

    print("GC-content:", gc_percent)
