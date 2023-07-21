import numpy as np

AMINO_MAP = {'<pad>':24, '*': 23, 'A': 0, 'C': 4, 'B': 20,
             'E': 6, 'D': 3, 'G': 7, 'F': 13, 'I': 9, 'H': 8,
             'K': 11, 'M': 12, 'L': 10, 'N': 2, 'Q': 5, 'P': 14,
             'S': 15, 'R': 1, 'T': 16, 'W': 17, 'V': 19, 'Y': 18,
             'X': 22, 'Z': 21}
AMINO_MAP_REV = ['A','R','N','D','C','Q','E','G','H','I','L','K',
                 'M','F','P','S','T','W','Y','V','B','Z','X','*','@']

valid_amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K',
                 'M','F','P','S','T','W','Y','V','B','Z','X']

def filter_seqs(seqs):
    filtered_strings = []
    for seq in seqs:
        valid = True
        for char in seq:
            if char not in valid_amino_acids:
                valid = False
                break
        if valid:
            filtered_strings.append(seq)
    return filtered_strings

def pad(seqs, init_token=None, eos_token=None, pad_token="<pad>", stop_words=None, fix_length=None, pad_type="mid", truncate_first=False):
    try:
        stop_words = set(stop_words) if stop_words is not None else None
    except TypeError:
        raise ValueError("Stop words must be convertible to a set")

    if fix_length is None:
        max_len = max(len(x) for x in seqs)
    else:
        max_len = fix_length + (
            init_token, eos_token).count(None) - 2
    padded, lengths = [], []
    for x in seqs:
        if pad_type == 'front':
            padded.append(
                [pad_token] * max(0, max_len - len(x))
                + ([] if init_token is None else [init_token])
                + list(x[-max_len:] if truncate_first else x[:max_len])
                + ([] if eos_token is None else [eos_token]))
        elif pad_type == 'end':
            padded.append(
                ([] if init_token is None else [init_token])
                + list(x[-max_len:] if truncate_first else x[:max_len])
                + ([] if eos_token is None else [eos_token])
                + [pad_token] * max(0, max_len - len(x)))
        elif pad_type == 'mid':
            i_gap = np.int32(np.ceil(min(len(x), max_len) / 2))
            i_gap_rev = min(len(x), max_len) - i_gap
            padded.append(
                ([] if init_token is None else [init_token])
                + list(x[:i_gap])
                + [pad_token] * max(0, max_len - len(x))
                + list(x[-i_gap_rev:])
                + ([] if eos_token is None else [eos_token]))
        else:
            raise ValueError('pad_type should be "front", "mid", or "end"')

        lengths.append(len(padded[-1]) - max(0, max_len - len(x)))

    return padded

def numerialize(seqs):

    numerialize_seqs = []

    for seq in seqs:
        encoded_text = [AMINO_MAP[char] for char in seq]
        numerialize_seqs.append(encoded_text)

    return numerialize_seqs
