def avg_phred_pos10(filename):
 
    total = 0
    count = 0

    with open(filename, "r", encoding="utf-8") as f:

        while True:

            header = f.readline()
            if not header:
                break
            #парсим
            sequence = f.readline().strip()
            plus = f.readline()
            quality = f.readline().strip()

            # отсеиваем прочтения где нет 10-ой позиции (короче 10)
            if len(quality) >= 10:

                q_char = quality[9]
                # используем нашу формулу
                q = ord(q_char) - 33
                total += q
                count += 1
    avg_quality = round(total_quality / count)

    return avg_quality


if __name__ == "__main__":

    fname = "reads.fastq"
    res = avg_phred_pos10(fname)
    print(res)
