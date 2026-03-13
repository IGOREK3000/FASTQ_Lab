# если в посл-ти есиь N, то её качество считаем 0.
def get_qualities_with_zero_ns(seq, qual):
    q = []
    for el, q_char in zip(seq, qual):
        if el == "N":
            q.append(0)
        else:
            q.append(ord(q_char) - 33)
    return q

# None, если read полностью удаляется и последовательность с качеством, если остается (как в исходном коде)
def sliding_window_trim(seq, qual, window_length=5, required_quality=30):
    quals = get_qualities_with_zero_ns(seq, qual)
    total_required_quality = required_quality * window_length

    # если read короче окна, то удаляется
    if len(quals) < window_length:
        return None
    # первое окно 
    total = sum(quals[:window_length])
    if total < total_required_quality:
        return None
    #сколько первых символов последовательности мы пока планируем оставит (изначально весь read)
    length_to_keep = len(quals)
    # скользящее окно
    for i in range(len(quals) - window_length):
         # убираем левый символ и добавляем новый справа
        total = total - quals[i] + quals[i + window_length]

        if total < total_required_quality:
            # на веремя ставим границу в правый край этого окна
            length_to_keep = i + window_length
            break

    # идни влево, пока последний сохраняемый символ не подходит
    i = length_to_keep
    last_base_quality = quals[i - 1]

    while last_base_quality < required_quality and i > 1:
        i -= 1
        last_base_quality = quals[i - 1]

    if i < 1:
        return None
    # если read стал короче (до i)
    if i < len(quals):
        return seq[:i], qual[:i]

    return seq, qual


def process_fastq(filename, window_length=5, required_quality=30, min_length=60):
    # сколько прочтений было полностью удалено триммингом
    fully_removed_by_quality = 0
    
    lengths_after_quality = []

    with open(filename, "r", encoding="utf-8") as f:
        while True:
            header = f.readline()
            if not header:
                break
            # парсим
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            trimmed = sliding_window_trim(
                seq,
                qual,
                window_length=window_length,
                required_quality=required_quality
            )
             # если read удалён полностью
            if trimmed is None:
                fully_removed_by_quality += 1
            else:
                trimmed_seq, trimmed_qual = trimmed
                lengths_after_quality.append(len(trimmed_seq))
    # минимальная средняя и максимальная длина после quality trimming
    min_len = min(lengths_after_quality)
    avg_len = round(sum(lengths_after_quality) / len(lengths_after_quality))
    max_len = max(lengths_after_quality)
    # фильтруем по длине уже после тримминга (второй фильтр по длине 60)
    remaining_after_two_filters = sum(
        1 for length in lengths_after_quality if length >= min_length
    )

    return (
        fully_removed_by_quality,
        min_len,
        avg_len,
        max_len,
        remaining_after_two_filters
    )


if __name__ == "__main__":
    fname = "reads.fastq"

    fully_removed, min_len, avg_len, max_len, remaining = process_fastq(
        fname,
        window_length=5,
        required_quality=30,
        min_length=60
    )

    print(fully_removed)
    print(min_len)
    print(avg_len)
    print(max_len)
    print(remaining)
